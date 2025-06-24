import os
import csv
import argparse
import yaml
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pymol import cmd
import requests
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform

def fetch_pdbs_by_uniprot(uniprot_id, output_folder):
    print(f"\n Поиск PDB структур для UniProt ID: {uniprot_id}")
    try:
        url = f"https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{uniprot_id}"
        response = requests.get(url)
        if response.status_code != 200:
            print("❌ Не удалось получить данные с PDBe API.")
            return []

        pdb_list = response.json().get(uniprot_id, [])
        pdb_ids = sorted(set(entry["pdb_id"].upper() for entry in pdb_list))  # ✅ Убираем дубликаты

        print(f"✅ Найдено {len(pdb_ids)} уникальных PDB структур: {', '.join(pdb_ids)}")

        downloaded = []
        for pdb_id in pdb_ids:
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            dest_path = os.path.join(output_folder, f"{pdb_id}.pdb")

            if not os.path.exists(dest_path):
                r = requests.get(url)
                if r.status_code == 200:
                    with open(dest_path, 'w') as f:
                        f.write(r.text)
                    print(f"⬇️  Загружено: {pdb_id}")
                    downloaded.append(pdb_id)
                else:
                    print(f"⚠️  Не удалось загрузить: {pdb_id}")
            else:
                print(f" Уже существует: {pdb_id}")
                downloaded.append(pdb_id)

        return downloaded

    except Exception as e:
        print(f"❌ Ошибка при загрузке: {e}")
        return []

def fetch_uniprot_by_pdb(pdb_id):
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"
    response = requests.get(url)
    if response.status_code != 200:
        return None
    data = response.json()
    mappings = data.get(pdb_id.lower(), {}).get("UniProt", {})
    if mappings:
        return list(mappings.keys())[0]
    return None

def compute_rmsd(sel1, sel2, method="super"):
    n1 = cmd.count_atoms(sel1)
    n2 = cmd.count_atoms(sel2)

    if n1 == 0 or n2 == 0:
        print(f"⚠️ Пропуск RMSD между {sel1.replace('pocket_', '')} | {sel2.replace('pocket_', '')} — пустая селекция.")
        return None, f"Пустая селекция: {n1} vs {n2}"
    
    allow_mismatch = int(n1 * 0.5)
    if abs(n1 - n2) > allow_mismatch:
        print(f"⚠️ Пропуск RMSD между {sel1} и {sel2} — слишком разное количество атомов ({n1} vs {n2}).")
        return None, f"Слишком разное количество атомов: {n1} vs {n2}"

    try:
        if method in ["rms", "rms_cur"] and n1 != n2:
            msg = f"⚠️ Метод {method} требует равное количество атомов: {n1} vs {n2}"
            print(msg)
            return None, msg

        if method == "align":
            rmsd = cmd.align(sel1, sel2)[0]
        elif method == "super":
            rmsd = cmd.super(sel1, sel2)[0]
        elif method == "rms":
            rmsd = cmd.rms(sel1, sel2)
        elif method == "rms_cur":
            rmsd = cmd.rms_cur(sel1, sel2)
        else:
            raise ValueError(f"Неизвестный метод RMSD: {method}")

        return rmsd, f"Метод: {method}, атомов: {n1} vs {n2}"

    except Exception as e:
        print(f"❌ Ошибка при расчёте RMSD методом {method} между {sel1} и {sel2}: {e}")
        return None, f"{e} | Метод: {method}, атомов: {n1} vs {n2}"

def print_directory_contents_pretty(path, columns=3):
    items = sorted(os.listdir(path))
    print(f"\n Рабочая директория: {path}\n")

    # Добавляем иконки и сортируем
    pretty_items = []
    for name in items:
        full_path = os.path.join(path, name)
        if os.path.isdir(full_path):
            pretty_items.append(f"[#] {name}")
        else:
            pretty_items.append(f"[*] {name}")

    # Печатаем в 3 колонки
    print(" Содержимое директории:\n")
    for i in range(0, len(pretty_items), columns):
        row = pretty_items[i:i+columns]
        print("   ".join(f"{x:<40}" for x in row))  # <40 — выравнивание по ширине

def clean_structure(name, mode):
    if mode == "1":
        cmd.remove(f"{name} and solvent")
    elif mode == "2":
        metals = ["MG", "ZN", "CA", "MN", "FE", "NA", "K", "CU", "CO", "NI", 'CL']
        metal_selection = " or ".join([f"{name} and resn {ion}" for ion in metals])
        cmd.remove(f"{name} and solvent")
        cmd.remove(metal_selection)
    elif mode == "3":
        cmd.remove(f"{name} and not polymer")
    elif mode == "4":
        pass
    else:
        print("⚠️ Неверный режим очистки. Пропуск.")

def load_yaml_config(config_path):
    print(f"\nПопытка загрузки YAML-конфигурации из: {os.path.abspath(config_path)}")
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)  # Используем safe_load, а не load
            print(f"Конфигурация загружена:\n{config}\n")
            return config
    except Exception as e:
        print(f"Ошибка при загрузке YAML-файла: {e}\n")
        return None

# --- Список ионов металлов ---
ions = ["MG", "ZN", "CA", "MN", "FE", "NA", "K", "CU", "CO", "NI"]

def parse_args():
    parser = argparse.ArgumentParser(description="Скрипт для анализа PDB файлов с RMSD и кластеризацией.")
    parser.add_argument('--config', type=str, help="Путь к конфигурационному файлу (JSON).")
    
    parser.add_argument('--mode', type=str, choices=[1,2,3], help="Режим работы: 1 - локальная папка, 2 - UniProt ID, 3 - PDB ID")
    parser.add_argument('--folder_path', type=str, help="Путь к папке с PDB файлами")
    parser.add_argument('--uniprot_id', type=str, help="UniProt ID для загрузки PDB структур")
    parser.add_argument('--pdb_id', type=str, help="PDB ID для поиска UniProt ID")

    # parser.add_argument('--folder_path', type=str, help="Путь к папке с PDB файлами.")
    parser.add_argument('--do_preprocess', type=str, choices=['0', '1'], help="Провести предварительную очистку (1 - да, 0 - нет).")
    parser.add_argument('--clean_options', type=str, nargs='*', help="Опции очистки (1-6, 8 для сохранения).")
    parser.add_argument('--save_dir', type=str, help="Папка для сохранения обработанных структур.")
    parser.add_argument('--ref_pocket', type=int, help="Индекс референсной структуры для кармана (1-based).")
    parser.add_argument('--ref_align', type=int, help="Индекс референсной структуры для выравнивания (1-based).")
    parser.add_argument('--pocket_method', type=str, choices=['1', '2', '3', '4'], help="Метод задания кармана.")
    parser.add_argument('--ligand_resi', type=str, help="Номер остатка (например, '40').")
    parser.add_argument('--ligand_chain', type=str, help="chain ID (например, 'A').")
    parser.add_argument('--radius', type=float, help="Радиус для кармана (в Å).")
    parser.add_argument('--resi_chain', type=str, help="Список остатков для метода 2 (например, 'PRO 46 A, ASN 61 A').")
    parser.add_argument('--chain_id', type=str, help="Идентификатор цепи для метода 3.")
    parser.add_argument('--comparison_mode', type=str, choices=['1', '2'], help="Режим сравнения (1 - all vs all, 2 - all vs ref).")
    parser.add_argument('--rmsd_method', type=str, choices=['1', '2', '3', '4'], help="Метод расчёта RMSD.")
    parser.add_argument('--linkage_method', type=str, choices=['1', '2', '3', '4', '5', '6', '7'], help="Метод кластеризации.")
    parser.add_argument('--cl_choice', type=str, choices=['1', '2', '3'], help="Способ кластеризации.")
    parser.add_argument('--num_clusters', type=int, help="Количество кластеров (для cl_choice 1).")
    parser.add_argument('--distance_threshold', type=float, help="Порог расстояния (для cl_choice 2).")
    return parser.parse_args()

# Очистка среды PyMOL
cmd.reinitialize()

print("\n----------------------------------------------------------------------------------------------------")
print("                                        Загрузка Аргументов                                         ")
print("----------------------------------------------------------------------------------------------------")

args = parse_args()
config = None
if args.config:
    config = load_yaml_config(args.config)
    if not config:
        print("Не удалось загрузить конфигурацию. Переход в интерактивный режим.")
        config = {}
else:
    config = {}

# # Заполняем недостающие аргументы значениями из конфигурационного файла
# for key, value in config.items():
#     if getattr(args, key, None) is None:
#         setattr(args, key, value)

print(f"\n Аргументы - args: {args}")
print(f"\n Аргументы из Конфигурации: {config}")

print("\n----------------------------------------------------------------------------------------------------")
print("                                        Загрузка PDB файлов                                         ")
print("----------------------------------------------------------------------------------------------------")

mode = args.mode or config.get('mode')
folder_path = args.folder_path or config.get('folder_path')
uniprot_id = args.uniprot_id or config.get('uniprot_id')
pdb_id = args.pdb_id or config.get('pdb_id')

# Показываем содержимое текущей директории
print_directory_contents_pretty(os.getcwd())

if not mode:
    # интерактивный ввод
    print("\nВыбери режим работы:")
    print("1 — Использовать локальную папку с PDB файлами")
    print("2 — Скачать все PDB структуры по UniProt ID")
    print("3 — Скачать все PDB структуры по PDB ID → UniProt ID")
    while True:
        mode = input("Введи номер режима (1/2/3): ").strip()
        if mode not in {"1", "2", "3"}:
            print("⚠️ Неверный ввод. Попробуй снова.")
        else:
            mode = int(mode)
            break

mode = str(mode)

if mode == "1":
    if not folder_path:
        print_directory_contents_pretty(os.getcwd())
        folder_path = input("\nВведи путь к папке с PDB файлами (по умолчанию текущая): ").strip()
        if not folder_path:
            folder_path = os.getcwd()

elif mode == "2":
    if not uniprot_id:
        uniprot_id = input("🔹 Введи UniProt ID: ").strip()
    folder_path = os.path.join(os.getcwd(), f"{uniprot_id}_pdbs")
    os.makedirs(folder_path, exist_ok=True)
    dwn=fetch_pdbs_by_uniprot(uniprot_id, folder_path)
    if dwn == []:
        print("❌ Не удалось найти UniProt ID")
        os.rmdir(folder_path)
        exit()
    print(f"✅ UniProt ID: {uniprot_id}")


elif mode == "3":
    if not pdb_id:
        pdb_id = input("🔹 Введи PDB ID: ").strip().lower()
    uniprot_id = fetch_uniprot_by_pdb(pdb_id)
    if not uniprot_id:
        print("❌ Не удалось найти UniProt ID по данному PDB ID.")
        exit()
    print(f"✅ UniProt ID: {uniprot_id}")
    folder_path = os.path.join(os.getcwd(), f"{uniprot_id}_pdbs")
    os.makedirs(folder_path, exist_ok=True)
    fetch_pdbs_by_uniprot(uniprot_id, folder_path)

# Проверяем, что папка существует и читаем PDB файлы
if not folder_path or not os.path.exists(folder_path):
    print(f"❌ Папка {folder_path} не найдена или не указана.")
    exit()

files = [f for f in os.listdir(folder_path) if f.endswith('.pdb')]
files.sort()

print(f"\n✅ Файлы успешно загружены из директории: {folder_path}")

print("\n----------------------------------------------------------------------------------------------------")
print("                                       Предварительная очистка                                      ")
print("----------------------------------------------------------------------------------------------------")

# --- Запрос на предварительную очистку ---
do_preprocess = args.do_preprocess or config.get('do_preprocess')
if do_preprocess is not None:
    do_preprocess = str(do_preprocess)
if not do_preprocess:
    do_preprocess = input("\nПровести предварительную очистку структур? (1 - да, 0 - нет): ").strip()
    if do_preprocess == '1':
        print("\nПожалуйста, выберите подходящие опции предварительной обработки структур:")
        print("1 — Удалить воду (solvent)")
        print("2 — Удалить ионы металлов (metal ions)")
        print("3 — Удалить всё, кроме белка (оставить только полимерную цепь)")
        print("4 — Удалить альтернативные конформации (altloc)")
        print("5 — Удалить атомы водорода (H)")
        print("6 — Удалить анизотропные параметры (ANISOU)")
        print("7 — Не очищать / Завершить выбор")
        print("8 - Сохранить обработанные структуры в определённой папке")
        
        save_dir = None
        to_remove = set()
        
        while True:
            choice = input("\nВведите номер опции (1–7), чтобы выбрать действие, или 7 — для завершения: ").strip()
            if choice == "7":
                break
            elif choice in ["1", "2", "3", "4", "5", "6"]:
                to_remove.add(choice)
                print(f"Опция {choice} будет применена.")
            elif choice == "8":
                save_dir = input("Введите путь к папке для сохранения обработанных структур: ").strip()
                if not os.path.exists(save_dir):
                    try:
                        os.makedirs(save_dir)
                        print(f"Папка не найдена. Создаю новую директорию: {os.path.abspath(save_dir)}")
                    except Exception as e:
                        print(f"Не удалось создать папку. Ошибка: {e}")
                        save_dir = None
                else:
                    print(f"Структуры будут сохранены в существующую папку: {os.path.abspath(save_dir)}")
            else:
                print("Некорректный ввод. Пожалуйста, введите число от 1 до 7")

        print("\nВыбранные опции будут применены: " + (", ".join(sorted(to_remove)) if to_remove else "ничего"))
else:
    save_dir = args.save_dir or config.get('save_dir')
    print(f"\nСтруктуры будут сохранены в существующую папку: {os.path.abspath(save_dir)}")
    to_remove = args.clean_options or config.get('clean_options') or []
    to_remove = set(to_remove)
    print("Выбранные опции будут применены: " + (", ".join(sorted(to_remove)) if to_remove else "ничего"))


# --- Загрузка и обработка всех структур ---
names = []
for f in files:
    name = os.path.splitext(f)[0]
    names.append(name)
    cmd.load(os.path.join(folder_path, f), name)


    # --- Применение выбранных опций очистки ---
    if do_preprocess == '1':
        if "1" in to_remove:
            cmd.remove(f"{name} and solvent")
        if "2" in to_remove:
            for ion in ions:
                cmd.remove(f"{name} and resn {ion}")
        if "3" in to_remove:
            cmd.remove(f"{name} and not polymer")
        if "4" in to_remove:
            cmd.remove(f"{name} and not alt '' and not alt A")
        if "5" in to_remove:
            cmd.remove(f"{name} and elem H")
        if "6" in to_remove:
            cmd.remove(f"{name} and anisou")

    # Сохраняем обработанную структуру, если указан путь
    if save_dir is not None:
        save_path = os.path.join(save_dir, f"{name}_processed.pdb")
        cmd.save(save_path, name)
        print(f"Сохранено: {save_path}")

print("\n----------------------------------------------------------------------------------------------------")
print("                                      Выбор референсных структур                                    ")
print("----------------------------------------------------------------------------------------------------")

# --- Выбор референсных структур ---
ref_name = args.ref_pocket or config.get('ref_pocket')

if not ref_name:
    print("\nДоступные PDB файлы:")
    for i, f in enumerate(files, 1):
        print(f"{i}. {f}")

    ref_index = int(input("\nВыберите референсную структуру для определения кармана: ")) - 1
    ref_name = names[ref_index]
else:
    print(f"\nРеференснуая структура для определения кармана: {ref_name}")

ref_align_name = args.ref_align or config.get('ref_align')

if not ref_align_name:
    print("\nДоступные PDB файлы:")
    for i, f in enumerate(files, 1):
        print(f"{i}. {f}")
    ref_align_index = int(input("\nВыберите референсную структуру для выравнивания: ")) - 1
    ref_align_name = names[ref_align_index]
else:
    print(f"\nРеференсная структура для выравнивания: {ref_align_name}")


# список всех HET-групп в выбранной структуре ref_index
het_atoms = cmd.get_model(f"{ref_name} and not polymer and not solvent").atom
het_residues = sorted({(a.resn, a.resi, a.chain) for a in het_atoms}, key=lambda x: (int(x[1]), x[2]))

print("\n----------------------------------------------------------------------------------------------------")
print("                                   Выбор метода задания кармана                                     ")
print("----------------------------------------------------------------------------------------------------")


# --- Выбор метода задания кармана ---
pocket_method = args.pocket_method or config.get('pocket_method')
if not pocket_method:
    # Показать список всех HET-групп 
    if het_residues:
        print(f" \nОбнаружены HET-группы в структуре {ref_name}:")
        for resn, resi, chain in het_residues:
            print(f"  - {resn:>3} {resi:>4} {chain}")
    else:
        print(f"\nВ структуре {ref_name} не обнаружено HET-групп.")
    
    print(f"\nВыбери способ задания кармана:")
    print("1 - Радиальное ( ввести остаток ID и радиус в Å)")
    print("2 - Вручную (ввести список остатков)")
    print("3 - По всей цепи (указать chain ID)")
    print("4 - По всем HET-группам во всех структурах (в заданном радиусе Å)")
    pocket_method = input("\nВведи способ (1, 2, 3 или 4): ").strip()

    if pocket_method == "1":
        while True:
            ligand_input = input("\nВведи номер остатка и chain ID через пробел, например '40 A': ").strip()
            parts = ligand_input.split()
            if len(parts) != 2:
                print("Ошибка: нужно ввести номер остатка и chain ID через пробел, например '40 A'. Попробуй ещё раз.")
                continue
            ligand_resi, ligand_chain = parts
            if not ligand_resi.isdigit():
                print("Ошибка: номер остатка должен быть числом. Попробуй ещё раз.")
                continue
            break
        try:
            radius = float(input("\nВведи радиус для кармана (в Å, по умолчанию 7): ").strip())
        except ValueError:
            radius = 7.0
    elif pocket_method == "2":
        resi_chain_input = input("\nВведи список остатков в формате (resn resi chain), например: PRO 46 A, ASN 61 A: ").strip()
    elif pocket_method == "3":
        chain_id = input("\nВведи идентификатор цепи (chain ID), например: A: ").strip().upper()
    elif pocket_method == "4":
        try:
            radius = float(input("\nВведи радиус для карманов вокруг всех HET-групп (в Å, по умолчанию 7): ").strip())
        except ValueError:
            radius = 7.0
else:
    print(f"\nМетод задания кармана - {pocket_method}")
    ligand_resi = args.ligand_resi or config.get('ligand_resi')
    print(f"ligand_resi - {ligand_resi}")
    ligand_chain = args.ligand_chain or config.get('ligand_chain')
    print(f"ligand_chain - {ligand_chain}")
    radius = args.radius or config.get('radius', 8.0)
    print(f"radius - {radius}")
    resi_chain_input = args.resi_chain or config.get('resi_chain')
    print(f"resi_chain - {resi_chain_input}")
    chain_id = args.chain_id or config.get('chain_id')
    print(f"chain_id - {chain_id}")
    
pocket_method =str(pocket_method)   
print("\n----------------------------------------------------------------------------------------------------")
print("                                        Режим сравнения                                             ")
print("----------------------------------------------------------------------------------------------------")

# --- Режим сравнения ---
comparison_mode = args.comparison_mode or config.get('comparison_mode')
if not comparison_mode:
    print("\nВыбери режим сравнения:")
    print("1 - Все со всеми (all vs all)")
    print("2 - Все с референсом (all vs ref) ")
    comparison_mode = input("\nВведи режим (1 или 2): ").strip()
else:
    print(f"\ncomparison_mode - {comparison_mode}")
comparison_mode = str(comparison_mode)

print("\n----------------------------------------------------------------------------------------------------")
print("                                          Метод RMSD                                                ")
print("----------------------------------------------------------------------------------------------------")

# --- Метод RMSD ---
rmsd_method_input = args.rmsd_method or config.get('rmsd_method')
if not rmsd_method_input:
    print("\nВыбери метод расчёта RMSD:")
    print("1 - align   (Строгая оценка, исключает несоответствия)")
    print("2 - super   (Гибкий подход, учитывает несовпадения)")
    print("3 - rms     (Точное RMSD, требует полное соответствие атомов)")
    print("4 - rms_cur (Быстрее, приближённое RMSD, тоже требует полное соответствие)")
    rmsd_method_input = input("\nВведи метод (1, 2, 3 или 4): ").strip()

rmsd_method_input =str(rmsd_method_input)
if rmsd_method_input == "1":
    rmsd_method = "align"
elif rmsd_method_input == "2":
    rmsd_method = "super"
elif rmsd_method_input == "3":
    rmsd_method = "rms"
elif rmsd_method_input == "4":
    rmsd_method = "rms_cur"
else:
    raise ValueError("Неверный метод расчёта RMSD.")

print(f"\nrmsd_method - {rmsd_method}")

print("\n----------------------------------------------------------------------------------------------------")
print("                                    Метод кластеризации                                             ")
print("----------------------------------------------------------------------------------------------------")

# --- Метод кластеризации ---

linkage_method_choice = args.linkage_method or config.get('linkage_method')
if not linkage_method_choice:
    print("\nВыбери метод иерархической кластеризации:")
    print("1 - ward (минимизация внутрикластерной дисперсии, требует евклидово расстояние)")
    print("2 - single (минимальное расстояние между кластерами)")
    print("3 - complete (максимальное расстояние между кластерами)")
    print("4 - average (среднее расстояние между кластерами (UPGMA))")
    print("5 - centroid (расстояние между центрами масс кластеров)")
    print("6 - median (медианное расстояние между кластерами)")
    print("7 - weighted (взвешенное среднее расстояние (WPGMA))")

    help_table = """
    | Метод        | Компактные кластеры | Вытянутые кластеры | Чувствительность к выбросам |
    | ------------ | ------------------- | ------------------ | --------------------------- |
    | 1. ward      | ✅ отлично          | ❌ плохо          | ⚠️ средняя                    |
    | 2. single    | ❌ плохо            | ✅ отлично        | ⚠️ высокая                    |
    | 3. complete  | ✅ хорошо           | ❌ плохо          | ✅ устойчив                  |
    | 4. average   | ✅ универсально     | ✅ неплохо        | ⚠️ средняя                    |
    | 5. centroid  | ⚠️ нестабильно       | ⚠️ нестабильно     | ⚠️ нестабильно                |
    | 6. median    | ⚠️ нестабильно       | ⚠️ нестабильно     | ⚠️ нестабильно                |
    | 7. weighted  | ✅ хорошо           | ✅ хорошо         | ⚠️ средняя                    |
    """

    print(help_table)

    linkage_method_choice = input("\nВведи номер метода (1-7): ").strip()

method_map = {
    "1": "ward",
    "2": "single",
    "3": "complete",
    "4": "average",
    "5": "centroid",
    "6": "median",
    "7": "weighted"
}

if linkage_method_choice in method_map:
    linkage_method = method_map[linkage_method_choice]
    print(f"\nВыбран метод кластеризации: {linkage_method}")
else:
    linkage_method = "ward"  # метод по умолчанию
    print(f"\nНеверный ввод, выбран метод по умолчанию: {linkage_method}")

print("\n----------------------------------------------------------------------------------------------------")
print("                                      Способ кластеризации                                          ")
print("----------------------------------------------------------------------------------------------------")

# --- Способ кластеризации ---
cl_choice = args.cl_choice if args.cl_choice else config.get('cl_choice')
if not cl_choice:
    print("\nВыберите способ кластеризации:")
    print("1 — Задать количество кластеров (maxclust)")
    print("2 — Задать порог расстояния (distance)")
    print("3 — Автоматический выбор числа кластеров (альтернативный вариант)")

    cl_choice = input("\nВаш выбор (1, 2 или 3): ").strip()

    if cl_choice == '1':
        try:
            num_clusters = int(input("\nВведите желаемое количество кластеров: ").strip())
        except ValueError:
            print("Некорректное число кластеров.")
    elif cl_choice == '2':
        try:
            distance_threshold = float(input("\nВведите порог расстояния (например, 3.0): ").strip())
        except ValueError:
            print("Некорректный порог расстояния.")
    elif cl_choice == '3':
        # Пример альтернативного варианта: автоматический выбор числа кластеров по "максимальному скачку"
        # Здесь простой вариант — ищем максимальный скачок в расстояниях linkage и делим на это
        print(f"Автоматически будет выбран порог расстояния ! ")
    else:
        print("Некорректный выбор.")
else:
    print(f"\nВыбран Способ кластеризации: {cl_choice}")
    num_clusters = args.num_clusters or config.get('num_clusters')
    print(f"num_clusters: {num_clusters}")
    distance_threshold = args.distance_threshold or config.get('distance_threshold')
    print(f"num_clusters: {distance_threshold}")

print("\n----------------------------------------------------------------------------------------------------\n")

print("\n" + "-" * 100)
print( " " * 40 + "Запуск скрипта !" +  " " * 40)
print("-" * 100 + "\n")


# --- Получение кармана ---

if pocket_method == "1":

    lig_sel = f"{ref_name} and resi {ligand_resi} and chain {ligand_chain}"
    cmd.select("ligand_ref", lig_sel)
    cmd.select("pocket_ref", f"(byres (ligand_ref around {radius})) and {ref_name} and polymer")
    resid_list = list({(atom.resn, atom.resi, atom.chain) for atom in cmd.get_model("pocket_ref").atom})
elif pocket_method == "2":
    resid_list = []
    for res in resi_chain_input.split(","):
        res = res.strip()
        if not res:
            continue
        try:
            resn, resi, chain = res.split()
            resid_list.append((resn, resi, chain))
        except ValueError:
            print(f"⚠️ Пропуск некорректного формата: {res}")
    if not resid_list:
        raise ValueError("Не указаны валидные остатки.")

elif pocket_method == "3":
    chain_sel = f"{ref_name} and chain {chain_id} and polymer"
    resid_list = list({(atom.resn, atom.resi, atom.chain) for atom in cmd.get_model(chain_sel).atom})

elif pocket_method == "4":

    all_het_selections = []

    for name in names:
        het_atoms = cmd.get_model(f"{name} and not polymer and not solvent").atom
        unique_residues = set((a.resn, a.resi, a.chain) for a in het_atoms)

        for resn, resi, chain in unique_residues:
            lig_sel = f"{name} and resn {resn} and resi {resi} and chain {chain}"
            pocket_sel = f"(byres ({lig_sel} around {radius})) and {name} and polymer"
            # printov heto kareliya sa eljokel grel infoi mej )
            # cmd.select("sel_pocket_sel", pocket_sel)
            # mini_resid_list = list({(a.resn, a.resi, a.chain) for a in cmd.get_model("sel_pocket_sel").atom})
            all_het_selections.append(pocket_sel)

    cmd.select("pocket_ref", " or ".join(all_het_selections))

    resid_list = list({(a.resn, a.resi, a.chain) for a in cmd.get_model("pocket_ref").atom})
    resid_list.sort(key=lambda x: (x[2], int(x[1])))

    print(f"\nОбщий карман, сформированный по всем HET-группам (радиус {radius} Å):")
    for resn, resi, chain in resid_list:
        print(f"  - {resn:>3} {resi:>4} {chain}")


else:
    raise ValueError("Неверный метод задания кармана.")

resid_list.sort(key=lambda x: (int(x[1]), x[2]))
print(f"{ref_name} - Остатки кармана:", ', '.join([f"{resn} {resi} {chain}" for resn, resi, chain in resid_list]))
# Формирование PyMOL-селекции
pymol_sel_parts = [f"(resi {resi} and chain {chain})" for resn, resi, chain in resid_list]
pymol_sel = f"{ref_name} and polymer and (" + " or ".join(pymol_sel_parts) + ")"

# Команды для красивой визуализации
vis_code = f'''
# Выделение и окраска кармана
select pocket, {pymol_sel}
show sticks, pocket
color red, pocket

# Отображение Cα атомов как шариков
select pocket_CA, pocket and name CA
show spheres, pocket_CA
set sphere_scale, 0.3, pocket_CA
color tv_red, pocket_CA
set stick_radius, 0.15, selection=pocket


# Остальной белок — полупрозрачная поверхность
select rest_protein, {ref_name} and polymer and not pocket
show surface, rest_protein
set transparency, 0.6, rest_protein
color gray80, rest_protein
'''
print("\nКоманды PyMOL для красивой визуализации кармана:\n")
print(vis_code.strip())
print(" \n ")

# --- Создание селекций карманов и Cα ---
for name in names:
    parts = [f"(resi {resi} and chain {chain})" for resn, resi, chain in resid_list]
    pocket_sel = f"{name} and polymer and (" + " or ".join(parts) + ")"
    ca_sel = f"{pocket_sel} and name CA"
    cmd.select(f"pocket_{name}", pocket_sel)
    cmd.select(f"pocket_CA_{name}", ca_sel)

# --- Выравнивание всех по референсу ---
ref_ca = f"pocket_CA_{ref_align_name}"
os.makedirs(f"{folder_path}/aligned_output", exist_ok=True)
os.makedirs(f"{folder_path}/aligned_output/aligned_structures", exist_ok=True)

failed_align = []  # список для белков, которые не выровнялись

for name in names:
    if name == ref_align_name:
        continue
    moving_ca = f"pocket_CA_{name}"
    try:
        cmd.align(moving_ca, ref_ca)
        aligned_path = os.path.join(f"{folder_path}/aligned_output/aligned_structures", f"{name}_aligned_to_{ref_align_name}.pdb")
        cmd.save(aligned_path, name)
        print(f"✅ {name} сохранён как {name}_aligned_to_{ref_align_name}.pdb")
    except:
        print(f"❌ Не удалось выровнять {name} по {ref_align_name}")
        failed_align.append(name)  # Добавляем в список неудач


# --- Расчёт RMSD ---
rmsd_all = []
rmsd_ca = []

if comparison_mode == "1":
    # All vs all
    n = len(names)
    heatmap_matrix = np.full((n, n), np.nan)
    heatmap_labels = names

    for i in range(n):
        for j in range(n):
            if i == j:
                heatmap_matrix[i, j] = 0.0
                continue

            sel1_all = f"pocket_{names[i]}"
            sel2_all = f"pocket_{names[j]}"
            sel1_ca = f"pocket_CA_{names[i]}"
            sel2_ca = f"pocket_CA_{names[j]}"

            rms_all, reason_all = compute_rmsd(sel1_all, sel2_all, rmsd_method)
            rms_ca, reason_ca = compute_rmsd(sel1_ca, sel2_ca, rmsd_method)

            rms_all_s = f"{rms_all:.5f}" if rms_all is not None else ""
            rms_ca_s = f"{rms_ca:.5f}" if rms_ca is not None else ""

            # print(f"{names[i]} vs {names[j]} | RMSD (all): {rms_all_s} | RMSD (Cα): {rms_ca_s}")
            
            rmsd_all.append((names[i], names[j], rms_all_s, reason_all))
            rmsd_ca.append((names[i], names[j], rms_ca_s, reason_ca))

            if rms_all is not None:
                heatmap_matrix[i, j] = rms_all
                heatmap_matrix[j, i] = rms_all  # симметрично

    # Построение тепловой карты
    plt.figure(figsize=(12, 10))

    # Вычисление среднего RMSD по строкам (игнорируя диагональ и nan)
    mean_rmsd = []
    for i in range(n):
        values = [heatmap_matrix[i, j] for j in range(n) if i != j and not np.isnan(heatmap_matrix[i, j])]
        avg = np.mean(values) if values else np.nan
        mean_rmsd.append(avg)

    # Добавление среднего как последнего столбца
    heatmap_matrix_with_mean = np.hstack((heatmap_matrix, np.array(mean_rmsd).reshape(-1, 1)))
    heatmap_labels_with_mean = heatmap_labels + ["Mean"]

    sns.heatmap(heatmap_matrix_with_mean, xticklabels=heatmap_labels_with_mean, yticklabels=heatmap_labels,
                cmap='coolwarm', annot=True, fmt=".2f", linewidths=0.5, cbar_kws={"label": "RMSD (Å)"})

    plt.title(f"RMSD (all atoms) Heatmap – Method: {rmsd_method}")
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()

    heatmap_path = os.path.join(folder_path, "aligned_output", "rmsd_heatmap.png")
    plt.savefig(heatmap_path, dpi=300)
    print(f"\nТепловая карта RMSD сохранена как: rmsd_heatmap.png")

    # Подготовка к Построению дендрограммы
    df = pd.DataFrame(heatmap_matrix, index=heatmap_labels, columns=heatmap_labels)


    # Удалим строки и столбцы, где больше 90% значений — NaN
    threshold = 0.9  # 90%

    row_nan_fraction = df.isna().mean(axis=1)
    col_nan_fraction = df.isna().mean(axis=0)

    df_filtered = df.loc[row_nan_fraction < threshold, df.columns[col_nan_fraction < threshold]]


    print("Максимальное значение RMSD:", df_filtered.values.max())
    print("Минимальное значение RMSD:", df_filtered.values.min())


    # valid_values = df_filtered.mask(df_filtered == 0).stack()  # исключим диагональ (0.0)
    # max_rmsd = valid_values.max()

    # df_filtered.replace([np.inf, -np.inf], np.nan, inplace=True)
    # df_filled = df.fillna(max_rmsd + 5.0)

    if df_filtered.shape[0] < 2:
        print("Недостаточно данных для кластеризации (менее двух структур).")
    else:
        condensed = squareform(df_filtered.values, checks=False)
        Z = linkage(condensed, method=f'{linkage_method}')

        clusters = None
        clustering_description = ""

        if cl_choice == '1':
            try:
                clusters = fcluster(Z, t=num_clusters, criterion='maxclust')
                clustering_description = f"{num_clusters}_clusters"
            except ValueError:
                print("Некорректное число кластеров.")
        elif cl_choice == '2':
            try:
                clusters = fcluster(Z, t=distance_threshold, criterion='distance')
                clustering_description = f"threshold_{distance_threshold}"
            except ValueError:
                print("Некорректный порог расстояния.")
        elif cl_choice == '3':
            # Пример альтернативного варианта: автоматический выбор числа кластеров по "максимальному скачку"
            # Здесь простой вариант — ищем максимальный скачок в расстояниях linkage и делим на это
            distances = Z[:, 2]
            max_gap_idx = distances.argmax()
            optimal_threshold = distances[max_gap_idx]

            clusters = fcluster(Z, t=optimal_threshold, criterion='distance')
            clustering_description = f"auto_threshold_{optimal_threshold:.2f}"

            print(f"Автоматически выбран порог расстояния: {optimal_threshold:.2f}")
        else:
            print("Некорректный выбор.")

        plt.figure(figsize=(10, 7))
        dendrogram(Z, labels=df_filtered.index.tolist())
        plt.title(f"Иерархическая кластеризация на основе RMSD ({linkage_method})")
        plt.tight_layout()
        dendro_path = os.path.join(folder_path, "aligned_output", f"RMSD_dendrogram_{linkage_method}.png")
        plt.savefig(dendro_path)
        plt.close()
        print(f"Дендрограмма сохранена в: RMSD_dendrogram_{linkage_method}.png")

        if clusters is not None:
            cluster_df = pd.DataFrame({
                'Protein': df_filtered.index,
                'Cluster': clusters
            })
            output_csv_path = os.path.join(folder_path, "aligned_output", f"cluster_assignments_{clustering_description}_{linkage_method}.csv")
            cluster_df.to_csv(output_csv_path, index=False)
            print(f"Результаты кластеризации сохранены в: cluster_assignments_{clustering_description}_{linkage_method}.csv")

elif comparison_mode == "2":
    # Все с референсом
    for name in names:
        if name == ref_align_name:
            continue
        sel1_all = f"pocket_{ref_align_name}"
        sel2_all = f"pocket_{name}"
        sel1_ca = f"pocket_CA_{ref_align_name}"
        sel2_ca = f"pocket_CA_{name}"

        rms_all, reason_all = compute_rmsd(sel1_all, sel2_all, rmsd_method)
        rms_ca, reason_ca = compute_rmsd(sel1_ca, sel2_ca, rmsd_method)
        
        # Форматирование RMSD для красивого вывода
        rms_all_s = f"{rms_all:.5f}" if rms_all is not None else ""
        rms_ca_s = f"{rms_ca:.5f}" if rms_ca is not None else ""

        rmsd_all.append((ref_align_name, name, rms_all_s, reason_all))
        rmsd_ca.append((ref_align_name, name, rms_ca_s, reason_ca))

        # print(f"{ref_align_name} vs {name} | RMSD (all): {rms_all_s} | RMSD (Cα): {rms_ca_s}")

else:
    raise ValueError("Неверный режим сравнения.")


# --- Сохранение RMSD ---
with open(f"{folder_path}/aligned_output/rmsd_all_atoms.csv", "w", newline='', encoding="utf-8-sig") as f:
    writer = csv.writer(f)
    writer.writerow(["Protein1", "Protein2", "RMSD_all_atoms", "Comment"])
    writer.writerows(rmsd_all)

with open(f"{folder_path}/aligned_output/rmsd_calpha.csv", "w", newline='', encoding="utf-8-sig") as f:
    writer = csv.writer(f)
    writer.writerow(["Protein1", "Protein2", "RMSD_Calpha", "Comment"])
    writer.writerows(rmsd_ca)

# --- Сохранение информации о карманах ---
with open(f"{folder_path}/aligned_output/info.txt", "w", encoding="utf-8-sig") as info_file:
    info_file.write(f"Структура для выравнивания: {ref_align_name}\n")
    if pocket_method == "1":
        info_file.write(f"Метод определения кармана: Радиальное (радиус {radius} Å от остатка {ligand_resi} {ligand_chain})\n")
    elif pocket_method == "2":
        info_file.write(f"Метод определения кармана: вручную (введённый список остатков)\n")
    elif pocket_method == "3":
        info_file.write(f"Метод определения кармана: по цепи {chain_id}\n")
    elif pocket_method == "4":
        info_file.write(f"Метод определения кармана: По всем HET-группам во всех структурах (в заданном радиусе)")        
    
    info_file.write(f"метод расчёта RMSD: {rmsd_method}\n")

    if failed_align:
        info_file.write("\nБелки, которые не удалось выровнять:\n")
        for f_name in failed_align:
            info_file.write(f"- {f_name}\n")
    else:
        info_file.write("\n\nВсе белки успешно выровнены.\n")
    
    info_file.write(f"\n\nИнформация о карманах (resn resi chain):\n\n")
    for name in names:
        pocket_atoms = cmd.get_model(f"pocket_{name}").atom
        residues = {(atom.resi, atom.resn, atom.chain) for atom in pocket_atoms}
        residues = sorted(residues, key=lambda x: (int(x[0]), x[2]))
        residue_lines = [f"{resn} {resi} {chain}" for resi, resn, chain in residues]
        info_file.write(f"{name}:\n")
        info_file.write("  " + ', '.join(residue_lines) + "\n\n")

    # --- Сравнение карманов с референсом ---
    info_file.write("\nСравнение карманов с референсной структурой:\n\n")
    
    # Формируем словарь для референса: (resi, chain) → resn
    ref_residues = {
        (atom.resi, atom.chain): atom.resn
        for atom in cmd.get_model(f"pocket_{ref_align_name}").atom
    }
    
    for name in names:
        if name == ref_align_name:
            continue
        target_atoms = cmd.get_model(f"pocket_{name}").atom
        target_residues = {
            (atom.resi, atom.chain): atom.resn
            for atom in target_atoms
        }
    
        diffs = []
        for key in ref_residues:
            ref_resn = ref_residues[key]
            tgt_resn = target_residues.get(key)
            if tgt_resn is None:
                diffs.append(f"{key[0]} {key[1]}: {ref_resn} → отсутствует")
            elif tgt_resn != ref_resn:
                diffs.append(f"{key[0]} {key[1]}: {ref_resn} → {tgt_resn}")
    
        if diffs:
            info_file.write(f"{name}: отличается {len(diffs)} остатками:\n")
            for d in diffs:
                info_file.write(f"  {d}\n")
            info_file.write("\n")
        else:
            info_file.write(f"{name}: идентичен референсу по карману.\n\n")

# --- Построение гистограмм RMSD ---
try:
    import matplotlib.pyplot as plt
except ImportError:
    print("⚠️ Модуль matplotlib не найден. Установи его командой: pip install matplotlib")
    plt = None

try:
    import seaborn as sns
except ImportError:
    print("⚠️ Модуль seaborn не найден. Установи его командой: pip install seaborn")
    sns = None

# Устанавливаем стиль оформления
sns.set(style="whitegrid")

if plt and sns:
    def plot_rmsd_histogram(rmsd_values, title, filename):
        values = [float(val[2]) for val in rmsd_values if val[2]]
        if not values:
            print(f"⚠️ Нет данных для построения гистограммы: {filename}")
            return

        plt.figure(figsize=(10, 6))
        sns.histplot(values, bins=30, kde=True, color='cornflowerblue', edgecolor='black')

        # Вертикальная линия среднего значения RMSD
        mean_val = np.mean(values)
        plt.axvline(mean_val, color='red', linestyle='--', linewidth=2, label=f'Среднее = {mean_val:.2f} Å')
        plt.legend(fontsize=12)
        
        plt.title(title, fontsize=16, fontweight='bold')
        plt.xlabel("RMSD (Å)", fontsize=14)
        plt.ylabel("Частота", fontsize=14)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)

        # Добавим разметку на оси X по округлённым значениям RMSD
        min_val, max_val = min(values), max(values)
        plt.xticks([round(x, 2) for x in 
                    list(plt.xticks()[0]) 
                    if min_val <= x <= max_val])

        plt.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()
        plt.savefig(f"{folder_path}/aligned_output/{filename}", dpi=300)
        plt.close()
        print(f"Гистограмма сохранена: {filename}")

    plot_rmsd_histogram(rmsd_all, "Гистограмма RMSD (все атомы)", "rmsd_all_atoms_hist.png")
    plot_rmsd_histogram(rmsd_ca, "Гистограмма RMSD (Cα)", "rmsd_calpha_hist.png")


print("\n✅ Готово. Файлы сохранены:")
print("-----------------------------")
print(f"{folder_path}/aligned_output/")
print("- RMSD по всем атомам: rmsd_all_atoms.csv")
print("- RMSD по Cα: rmsd_calpha.csv")
print("- Остатки карманов (resn resi chain): info.txt")
print("- Выровненные структуры: папка aligned_structures/")
print(f"- Гистограмма RMSD (все атомы): rmsd_all_atoms_hist.png")
print(f"- Гистограмма RMSD (Cα): rmsd_calpha_hist.png")
if comparison_mode == "1": 
    print(f"- Тепловая карта: rmsd_heatmap.png")
    print(f"- Дендрограмма: RMSD_dendrogram_{linkage_method}.png ")
    print(f"- Результаты кластеризации: cluster_assignments_{clustering_description}_{linkage_method}.csv")
print("\n")