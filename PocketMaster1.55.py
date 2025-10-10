import os
import sys
import time
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
from kneed import KneeLocator
from scipy.cluster.hierarchy import linkage, leaves_list


ions = [
    # Щелочные металлы
    "LI", "NA", "K", "RB", "CS",
    # Щелочноземельные металлы
    "MG", "CA", "SR", "BA",
    # Переходные металлы
    "MN", "FE", "CO", "NI", "CU", "ZN",
    "CD", "HG", "PT", "AU", "AG",
    # Редкие переходные и лантаноиды
    "Y", "ZR", "MO", "RU", "RH", "PD", "W", "OS", "IR",
    "LA", "CE", "PR", "ND", "SM", "EU", "GD", "TB",
    "DY", "HO", "ER", "TM", "YB", "LU",
    # Актиноиды (редко в PDB, но встречаются)
    "U", "TH",
    # Галогениды
    "F", "CL", "BR", "I",
    # Другие анионы
    "SO4", "PO4", "NO3", "CO3", "MO4", "WO4", "SE4"
]


modified_residues = [
    "CSO",  # S-hydroxycysteine – окисленный цистеин (SG-OH)
    "MSE",  # Selenomethionine – селенометионин (используется при структурной биологии)
    "SEP",  # Phosphoserine – фосфорилированный серин
    "TPO",  # Phosphothreonine – фосфорилированный треонин
    "PTR",  # Phosphotyrosine – фосфорилированный тирозин

    "HIC",  # Hydroxyisoleucine или Hydroxylysine – гидроксилированный изолейцин/лизин
    "HYP",  # Hydroxyproline – гидроксипролин (часто в коллагене)
    "SEC",  # Selenocysteine – селенцистеин (21-я аминокислота)
    "PYL",  # Pyrrolysine – пирролизин (22-я аминокислота, встречается у архей)

    "CME",  # S-methylcysteine – метилированный цистеин
    "FME",  # N-formylmethionine – инициирующий метионин в прокариотах
    "CSS",  # Disulfide bond – формализованная запись дисульфидного моста
    "CSD",  # Dithiothreitol-modified cysteine – восстановленная форма модифицированного цистеина
    "CSX",  # Mixed disulfide – смешанные формы модифицированных цистеинов

    "KCX",  # Carboxyglutamate – γ-карбоксиглутамат (встречается в белках, связывающих кальций)
    "LLP",  # Lysylpyridinoline – сшитый лизин (в коллагене)
    "MLY",  # N-methyl-lysine – метилированный лизин
    "MLZ",  # Ди- или три-метилированный лизин
    "ALY",  # N-acetyl-lysine – ацетилированный лизин

    "PCA",  # Pyroglutamate – циклизованный глутамин/глутамат на N-конце
]

# Буферные компоненты
buffer_residues = [
    "TRS",  # Tris  
    "MES",  # 2-(N-Morpholino)ethanesulfonic acid  
    "HEP",  # HEPES  
    "ACT",  # Acetate  
    "MOP",  # MOPS  
    "PIP",  # PIPES  
    "CAP",  # CAPS  
    "BIS",  # Bis-Tris  
    "TES",  # TES  
    "ADA",  # ADA  
    "AEC",  # ACES  
]

# Криопротектанты
cryo_residues = [
    "GOL",  # Glycerol  
    "EDO",  # Ethylene glycol  
    "MPD",  # 2-Methyl-2,4-pentanediol  
    "DMS",  # DMSO  
    "BME",  # β-mercaptoethanol  
    "PGO",  # Propylene glycol  
    "EGD",  # Ethylene glycol dimer  
    "PRP",  # n-Propanol  
]

# Сульфаты / Фосфаты
sulfate_phosphate_residues = [
    "SO4",  # Sulfate  
    "PO4",  # Phosphate  
    "POM",  # Pyrophosphate  
    "P5P",  # Polyphosphate  
]

# Восстановители
reductant_residues = [
    "DTT",   # Dithiothreitol  
    "BME",   # β-mercaptoethanol  
    "TCEP",  # Tris(2-carboxyethyl)phosphine  
]


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
        if "_CA_" not in sel1:
            print(f"⚠️ Пропуск RMSD между {sel1.replace('pocket_', '')} | {sel2.replace('pocket_', '')} — пустая селекция.")
        return None, f"Пустая селекция: {n1} vs {n2}"
    
    # allow_mismatch = int(n1 * 0.9)
    # if abs(n1 - n2) > allow_mismatch:
    #     if "_CA_" not in sel1:
    #         print(f"⚠️ Пропуск RMSD между {sel1} и {sel2} — слишком разное количество атомов ({n1} vs {n2}).")
    #     return None, f"Слишком разное количество атомов: {n1} vs {n2}"

    try:
        if method in ["rms", "rms_cur"] and n1 != n2:
            msg = f"⚠️ Метод {method} требует равное количество атомов: {n1} vs {n2}"
            print(msg)
            return None, msg

        if method == "align":
            rmsd = cmd.align(sel1, sel2)[0]
        elif method == "cealign":
            result = cmd.cealign(sel1, sel2)
            rmsd = result['RMSD']
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
        ions_selection = " or ".join([f"{name} and resn {ion}" for ion in ions])
        cmd.remove(f"{name} and solvent")
        cmd.remove(ions_selection)
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
    parser.add_argument('--init_align', type=str, choices=['1', '2'], help="Метод init_align")
    parser.add_argument('--init_align_chain_id', type=str, help=" init_align_chain_id chain ID (например, 'A')")
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
    parser.add_argument('--rmsd_method', type=str, choices=['1', '2', '3', '4', '5'], help="Метод расчёта RMSD.")
    parser.add_argument('--linkage_method', type=str, choices=['1', '2', '3', '4', '5', '6', '7'], help="Метод кластеризации.")
    parser.add_argument('--cl_choice', type=str, choices=['1', '2', '3'], help="Способ кластеризации.")
    parser.add_argument('--num_clusters', type=int, help="Количество кластеров (для cl_choice 1).")
    parser.add_argument('--distance_threshold', type=float, help="Порог расстояния (для cl_choice 2).")
    return parser.parse_args()


def save_config_to_yaml(config_dict, path="run_config.yaml"):
    # Удаляем пары с None
    cleaned_config = {k: v for k, v in config_dict.items() if v is not None}
    try:
        with open(path, "w") as f:
            yaml.dump(cleaned_config, f, default_flow_style=False, allow_unicode=True, sort_keys=False)
        print(f"\n✅ Конфигурация сохранена в файл: {path}")
    except Exception as e:
        print(f"❌ Ошибка при сохранении конфигурации: {e}")


# def build_config_from_args(args):
#     return {
#         "mode": args.mode,
#         "folder_path": args.folder_path,
#         "uniprot_id": args.uniprot_id,
#         "pdb_id": args.pdb_id,
#         "save_dir": args.save_dir,
#         "do_preprocess": args.do_preprocess,
#         "clean_options": args.clean_options,
#         "ref_pocket": args.ref_pocket,
#         "ref_align": args.ref_align,
#         "pocket_method": args.pocket_method,
#         "ligand_resi": args.ligand_resi,
#         "ligand_chain": args.ligand_chain,
#         "radius": args.radius,
#         "resi_chain": args.resi_chain,
#         "chain_id": args.chain_id,
#         "comparison_mode": args.comparison_mode,
#         "rmsd_method": args.rmsd_method,
#         "linkage_method": args.linkage_method,
#         "cl_choice": args.cl_choice,
#         "num_clusters": args.num_clusters,
#         "distance_threshold": args.distance_threshold,
#     }

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

ref_align_name = None
init_align = None
init_align_chain_id = None
failed_align = None
resi_chain = None
chain_id = None
num_clusters = None
distance_threshold = None


print("\n----------------------------------------------------------------------------------------------------")
print("                                        Загрузка PDB файлов                                         ")
print("----------------------------------------------------------------------------------------------------")

mode = args.mode or config.get('mode')
folder_path = args.folder_path or config.get('folder_path')
uniprot_id = args.uniprot_id or config.get('uniprot_id')
pdb_id = args.pdb_id or config.get('pdb_id')

# # Показываем содержимое текущей директории
# print_directory_contents_pretty(os.getcwd())

if not mode:
    # интерактивный ввод
    print("\nВыбери режим работы:")
    print("1 — Использовать локальную папку с PDB файлами")
    print("2 — На основе UniProt ID скачать все соответствующие PDB структуры")
    print("3 — На основе PDB ID определить UniProt ID и скачать все соответствующие PDB структуры")
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
    while True:
        do_preprocess = input("\nПровести предварительную очистку структур? (1 - да, 0 - нет): ").strip()
        if do_preprocess not in {"1", "0"}:
            print("⚠️ Неверный ввод. Попробуй снова.")
        else:
            break

    if do_preprocess == '1':
        print("\nПожалуйста, выберите подходящие опции предварительной обработки структур:")
        print("1 — Удалить воду (solvent)")
        print("2 — Удалить ионы (ions)")
        print("3 — Удалить сульфаты и фосфаты (SO4, PO4, и др.)")
        print("4 — Удалить буферные компоненты (TRS, MES, HEP, и др.)")
        print("5 — Удалить криопротектанты (GOL, EDO, MPD, и др.)")
        print("6 — Удалить восстановители (DTT, BME, TCEP)")
        print("7 — Удалить всё воду, ионы, буферы, криопротектанты, фосфаты, восстановители")
        print("8 — Удалить модифицированные аминокислотные остатки (CSO, MSE, SEP, TPO, PTR и др.)")
        print("9 — Удалить всё, кроме белка (оставить только полимерную цепь)")
        print("10 — Удалить альтернативные конформации (altloc)")
        print("11 — Удалить анизотропные параметры (ANISOU)")
        print("12 — Удалить атомы водорода (H)")
        print("13 - Сохранить обработанные структуры в определённой папке")
        print("14 — Не очищать / Завершить выбор")

        save_dir = None
        to_remove = set()
        
        while True:
            choice = input("\nВведите номер опции (1–13), чтобы выбрать действие, или 14 — для завершения: ").strip()
            if choice == "14":
                break
            elif choice in ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"]:
                to_remove.add(choice)
                print(f"Опция {choice} будет применена.")
            elif choice == "13":
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
                print("Некорректный ввод. Пожалуйста, введите число от 1 до 14")

        print("\nВыбранные опции будут применены: " + (", ".join(sorted(to_remove)) if to_remove else "ничего"))
    else:
        save_dir = None
        to_remove = None

else:
    save_dir = args.save_dir or config.get('save_dir')
    if save_dir:
        print(f"\nСтруктуры после предварительной очистки будут сохранены в существующую папку: {os.path.abspath(save_dir)}")
    else:
        print("\n⚠️ Структуры после предварительной очистки не будут сохранены в отделную папку. Папка для сохранения (\"save_dir\") не указана.")

    to_remove = args.clean_options or config.get('clean_options') or []
    to_remove = set(to_remove)
    print("\nБудут применены следующие опции обработки структур: " + (", ".join(sorted(to_remove)) if to_remove else "--> ⚠️ Не была выбрана ни одна из опций предварительной обработки"))


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
            for resn in sulfate_phosphate_residues:
                cmd.remove(f"{name} and resn {resn}")
        if "4" in to_remove:
            for resn in buffer_residues:
                cmd.remove(f"{name} and resn {resn}")
        if "5" in to_remove:
            for resn in cryo_residues:
                cmd.remove(f"{name} and resn {resn}")
        if "6" in to_remove:
            for resn in reductant_residues:
                cmd.remove(f"{name} and resn {resn}")
        if "7" in to_remove:
            # Объединённый список "ненужных" лиганды
            all_small_residues = set(
                buffer_residues +
                cryo_residues +
                sulfate_phosphate_residues +
                reductant_residues)

            # Удалить воду
            cmd.remove(f"{name} and solvent")
            
            # Удалить ионы
            for ion in ions:
                cmd.remove(f"{name} and resn {ion}")
            
            # Удалить все другие "мелкие" молекулы
            for resn in all_small_residues:
                cmd.remove(f"{name} and resn {resn}")
        if "8" in to_remove:
            for resn in modified_residues:
                cmd.remove(f"{name} and resn {resn}")
        if "9" in to_remove:
            cmd.remove(f"{name} and not polymer")
        if "10" in to_remove:
            cmd.remove(f"{name} and not (alt '' + A)")
        if "11" in to_remove:
            ani_state = 0
        if "12" in to_remove:
            cmd.remove(f"{name} and elem H")

        # Сохраняем обработанную структуру, если указан путь
        if save_dir is not None:
            save_path = os.path.join(save_dir, f"{name}_processed.pdb")
            cmd.save(save_path, name, state=0)
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

    ref_index = int(input("\nВыберите референсную структуру для выравнивания, введя её порядковый номер: ")) - 1
    ref_name = names[ref_index]
else:
    print(f"\nРеференснуая структура для выравнивания: {ref_name}")




# список всех HET-групп в выбранной структуре ref_index
het_atoms = cmd.get_model(f"{ref_name} and not polymer and not solvent").atom
het_residues = sorted({(a.resn, a.resi, a.chain) for a in het_atoms}, key=lambda x: (int(x[1]), x[2]))

# print("\n----------------------------------------------------------------------------------------------------")
# print("                   Предварительное выравнивание всех структур относительно референса                ")
# print("----------------------------------------------------------------------------------------------------")

# # --- Выравнивания всех по референсу ---
# init_align = args.init_align or config.get('init_align')
# if not init_align:    
#     print(f"\nВыберите метод выравнивания всех структур с референсной")
#     print("1 - По цепи (указать chain ID)")
#     print("2 - Автоматический (Без указания цепи)")
#     init_align = input("\nВведи способ (1 или 2): ").strip()
    
#     if init_align == "1":
#         init_align_chain_id = input("\nВведи идентификатор цепи (chain ID), например: A ").strip().upper()
#     elif init_align == "2":
#         print("\nВыбран Автоматический")
# else:
#     init_align =str(init_align)   
#     if init_align == "1":
#         init_align_chain_id = args.init_align_chain_id or config.get('init_align_chain_id')
#         print(f"\nВыравнивания всех структур с референсом: 1 - По цепи {init_align_chain_id}")
#     elif init_align == "2":
#         print("\nВыбран Автоматический")

print("\n----------------------------------------------------------------------------------------------------")
print("                         Выбор способа определения участка для выравнивания                         ")
print("----------------------------------------------------------------------------------------------------")


# --- Выбор метода задания кармана ---
ligand_resi = None
ligand_chain = None
pocket_method = args.pocket_method or config.get('pocket_method')
if not pocket_method:
    # Показать список всех HET-групп 
    if het_residues:
        print(f" \nОбнаружены HET-группы в структуре {ref_name}:")
        for resn, resi, chain in het_residues:
            print(f"  - {resn:>3} {resi:>4} {chain}")
    else:
        print(f"\nВ структуре {ref_name} не обнаружено HET-групп.")
    
    print(f"\nВыбери способ определения участка для выравнивания: участок определяется ... ")
    print("1 – На реф. структуре по заданному ID остатка и радиусу (Å), затем он ищется и выравнивается во всех структурах")
    print("2 – После предварительного выравнивания всех структур между собой, для каждой структуры вокруг выбранного реф. остатка")
    print("3 – Для каждой структуры вокруг её HET-групп в пределах заданного радиуса (Å)")
    print("4 – По введённому пользователем списку остатков, затем ищется во всех структурах")
    print("5 – На реф. структуре по указанному идентификатору цепи, затем ищется и выравнивается во всех структурах")
    
    while True:
        pocket_method = input("\nВведи способ (1, 2, 3, 4 или 5): ").strip()
        if pocket_method == "1":
            while True:
                ligand_input = input("\nВведи номер остатка и chain ID через пробел, например '40 A': ").strip()
                parts = ligand_input.split()
                if len(parts) != 2:
                    print("Ошибка: нужно ввести номер остатка и chain ID через пробел, например '40 A': Попробуй ещё раз.")
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
            break
        elif pocket_method == "4":
            resi_chain = input("\nВведи список остатков в формате (resn resi chain), например: PRO 46 A, ASN 61 A: ").strip()
            break
        elif pocket_method == "5":
            chain_id = input("\nВведи идентификатор цепи (chain ID), например: A: ").strip().upper()
            break
        elif pocket_method == "3":
            try:
                radius = float(input("\nВведи радиус для карманов вокруг всех HET-групп (в Å, по умолчанию 7): ").strip())
            except ValueError:
                radius = 7.0
            chain_id = input("\nВведи идентификатор цепи (chain ID), например: A: ").strip().upper()
            break
        elif pocket_method == "2":

            print("\n----------------------------------------------------------------------------------------------------")
            print("                   Предварительное выравнивание всех структур относительно референса                ")
            print("----------------------------------------------------------------------------------------------------")

            print("\nДоступные PDB файлы:")
            for i, f in enumerate(files, 1):
                print(f"{i}. {f}")
            ref_align_index = int(input("\nВыберите референсную структуру для предварительного выравнивания, введя её порядковый номер: ")) - 1
            ref_align_name = names[ref_align_index]
        
            # --- Выравнивания всех по референсу ---    
            print(f"\nВыберите метод предварительного выравнивания всех структур с референсной")
            print("1 - По цепи (указать chain ID)")
            print("2 - По всей молекуле")
            init_align = input("\nВведи способ (1 или 2): ").strip()
            
            if init_align == "1":
                init_align_chain_id = input("\nВведи идентификатор цепи (chain ID), например: A ").strip().upper()
            elif init_align == "2":
                print("\nПо всей молекуле")

            while True:
                ligand_input = input("\nДля определения карманов введи номер референсого остатка и chain ID через пробел, например '40 A': ").strip()
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
            break
        else:
            print("Некорректный ввод. Пожалуйста, введите число от 1 до 5")


else:
    pocket_method =str(pocket_method)   
    if pocket_method == "1":
        print("1 – На реф. структуре по заданному ID остатка и радиусу (Å), затем он ищется и выравнивается во всех структурах")
        ligand_resi = args.ligand_resi or config.get('ligand_resi')
        print(f"ligand_resi - {ligand_resi}")
        ligand_chain = args.ligand_chain or config.get('ligand_chain')
        print(f"ligand_chain - {ligand_chain}")
        radius = args.radius or config.get('radius')
        print(f"radius - {radius}")
    elif pocket_method == "4":
        print("4 – По введённому пользователем списку остатков, затем ищется во всех структурах")
        resi_chain = args.resi_chain or config.get('resi_chain')
        print(f"resi_chain - {resi_chain}")
    elif pocket_method == "5":
        print("5 – На реф. структуре по указанному идентификатору цепи, затем ищется и выравнивается во всех структурах")
        chain_id = args.chain_id or config.get('chain_id')
        print(f"chain_id - {chain_id}")
    elif pocket_method == "3":
        print("3 – Для каждой структуры вокруг её HET-групп в пределах заданного радиуса (Å)")
        radius = args.radius or config.get('radius')
        chain_id = args.chain_id or config.get('chain_id')
        print(f"chain_id - {chain_id}")
        print(f"radius - {radius}")
    if pocket_method == "2":
        print("2 – После предварительного выравнивания всех структур между собой, для каждой структуры вокруг выбранного реф. остатка")
        ref_align_name = args.ref_align or config.get('ref_align')
        print(f"\nРеференсная структура для предварительного выравнивания: {ref_align_name}")
        init_align = args.init_align or config.get('init_align')
        init_align = str(init_align)   
        if init_align == "1":
            init_align_chain_id = args.init_align_chain_id or config.get('init_align_chain_id')
            print(f"\n Предварительное Выравнивание всех структур с референсом: 1 - По цепи {init_align_chain_id}")
        elif init_align == "2":
            print("\nПо всей молекуле")
        ligand_resi = args.ligand_resi or config.get('ligand_resi')
        print(f"ligand_resi - {ligand_resi}")
        ligand_chain = args.ligand_chain or config.get('ligand_chain')
        print(f"ligand_chain - {ligand_chain}")
        radius = args.radius or config.get('radius')
        print(f"radius - {radius}")
    
print("\n----------------------------------------------------------------------------------------------------")
print("                                        Режим сравнения                                             ")
print("----------------------------------------------------------------------------------------------------")

# --- Режим сравнения ---
comparison_mode = args.comparison_mode or config.get('comparison_mode')
if not comparison_mode:
    print("\nВыбери режим сравнения:")
    print("1 - Все со всеми (all vs all)")
    print("2 - Все с референсом (all vs ref)")
    comparison_mode = input("\nВведи режим (1 или 2): ").strip()
else:
    comparison_mode = str(comparison_mode)
    if comparison_mode =="1":
            print("\nВыбран Режим сравнения: 1 - Все со всеми (all vs all)")
    elif comparison_mode =="2":
            print("\nВыбран Режим сравнения: 2 - Все с референсом (all vs ref)")

print("\n----------------------------------------------------------------------------------------------------")
print("                                          Метод RMSD                                                ")
print("----------------------------------------------------------------------------------------------------")

# --- Метод RMSD ---
rmsd_method_input = args.rmsd_method or config.get('rmsd_method')
if not rmsd_method_input:
    print("\nВыбери метод расчёта RMSD:")
    print("1 - align   (Строгая оценка, исключает несоответствия)")
    print("2 - cealign (Структурное выравнивание на основе геометрии, эффективно даже при низком сходстве)")
    print("3 - super   (Гибкий подход, учитывает несовпадения)")
    print("4 - rms     (Точное RMSD, требует полное соответствие атомов)")
    print("5 - rms_cur (Быстрее, приближённое RMSD, тоже требует полное соответствие)")
    rmsd_method_input = input("\nВведи метод (1, 2, 3 или 4): ").strip()

rmsd_method_input =str(rmsd_method_input)
if rmsd_method_input == "1":
    rmsd_method = "align"
elif rmsd_method_input == "2":
    rmsd_method = "cealign"
elif rmsd_method_input == "3":
    rmsd_method = "super"
elif rmsd_method_input == "4":
    rmsd_method = "rms"
elif rmsd_method_input == "5":
    rmsd_method = "rms_cur"
else:
    raise ValueError("Неверный метод расчёта RMSD.")

print(f"\nВыбран метод расчёта RMSD: {rmsd_method_input} - {rmsd_method}")

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

linkage_method_choice = str(linkage_method_choice)
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
    print(f"\nВыбран метод кластеризации: {linkage_method_choice} - {linkage_method}")
else:
    linkage_method = "ward"  # метод по умолчанию
    print(f"\nНеверный ввод, выбран метод по умолчанию: 1 - {linkage_method}")

print("\n----------------------------------------------------------------------------------------------------")
print("                         Параметры кластеризации: количество или порог                              ")
print("----------------------------------------------------------------------------------------------------")

# --- Способ кластеризации ---
cl_choice = args.cl_choice if args.cl_choice else config.get('cl_choice')
if not cl_choice:
    print("\nВыберите способ кластеризации:")
    print("1 — Задать количество кластеров (maxclust)")
    print("2 — Задать порог расстояния (distance)")
    print("3 — Автоматический: Порог на основе 70% процента максимального расстояния слиянияя") 
    print("4 — Автоматический: Метод Изгиба или Локтя (Elbow Method)")

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
        print(f"Автоматически будет выбран порог расстояния ! ")
    elif cl_choice == '4':
        print(f"Автоматически будет выбран порог расстояния ! ")
    else:
        print("Некорректный выбор.")
else:
    cl_choice = str(cl_choice)
    if cl_choice == '1':
        print(f"\nВыбран Способ кластеризации: 1 — Задать количество кластеров (maxclust)")
        num_clusters = args.num_clusters or config.get('num_clusters')
        print(f"num_clusters: {num_clusters}")
    elif cl_choice == '2':
        print(f"\nВыбран Способ кластеризации: 2 — Задать порог расстояния (distance)")
        distance_threshold = args.distance_threshold or config.get('distance_threshold')
        print(f"distance_threshold: {distance_threshold}")
    elif cl_choice == "3":
        print(f"\nВыбран Способ кластеризации: 3 — Автоматический: Порог на основе 70% процента максимального расстояния слиянияя") 
    elif cl_choice == "4":
        print(f"\nВыбран Способ кластеризации: 4 — Автоматический: Метод Изгиба или Локтя (Elbow Method)") 

#===================================================================================================================================================================
print("\n" + "-" * 100)
print( " " * 40 + "Запуск скрипта !" +  " " * 40)
print("-" * 100 + "\n")
#===================================================================================================================================================================

for name in names:
    print(f"✅ {name} Запуск обработки ")
    time.sleep(0.02)

#--------------------------------------------------------------------
# text_list = [f"✅ {name} Запуск обработки" for name in names]

# for text in text_list:
#     for char in text:
#         sys.stdout.write(char)
#         sys.stdout.flush()
#         time.sleep(0.005)  # задержка между символами 
#     sys.stdout.write("\n")
#     sys.stdout.flush()
#     time.sleep(0.03)  # задержка между строками
#--------------------------------------------------------------------

# --- Выравнивание всех по референсу ---

os.makedirs(f"{folder_path}/output", exist_ok=True)
os.makedirs(f"{folder_path}/output/aligned_structures", exist_ok=True)

# --- Получение кармана ---

if pocket_method == "1":
    lig_sel = f"{ref_name} and resi {ligand_resi} and chain {ligand_chain}"
    cmd.select("ligand_ref", lig_sel)
    cmd.select("pocket_ref", f"(byres (ligand_ref around {radius})) and {ref_name} and polymer")
    resid_list = list({(atom.resn, atom.resi, atom.chain) for atom in cmd.get_model("pocket_ref").atom})
elif pocket_method == "4":
    resid_list = []
    for res in resi_chain.split(","):
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

elif pocket_method == "5":
    chain_sel = f"{ref_name} and chain {chain_id} and polymer"
    resid_list = list({(atom.resn, atom.resi, atom.chain) for atom in cmd.get_model(chain_sel).atom})

elif pocket_method == "3":

    for name in names:
        het_atoms = cmd.get_model(f"{name} and not polymer and not solvent and chain {chain_id}").atom
        unique_residues = set((a.resn, a.resi, a.chain) for a in het_atoms)

        # Фильтрация: оставляем только HET-группы, которые не входят в список модифицированных остатков
        filtered_residues = [
            (resn, resi, chain)
            for (resn, resi, chain) in unique_residues
            if resn not in modified_residues
        ]

        for resn, resi, chain in filtered_residues:
            lig_sel = f"{name} and resn {resn} and resi {resi} and chain {chain}"
            pocket_sel = f"(byres ({lig_sel} around {radius})) and {name} and polymer"
            ca_sel = f"{pocket_sel} and name CA"
            cmd.select(f"pocket_{name}", pocket_sel)
            cmd.select(f"pocket_CA_{name}", ca_sel)

    resid_list = list({(atom.resn, atom.resi, atom.chain) for atom in cmd.get_model(f"pocket_{name}").atom})

    #         all_het_selections.append(pocket_sel)

    # cmd.select("pocket_ref", " or ".join(all_het_selections))

    # resid_list = list({(a.resn, a.resi, a.chain) for a in cmd.get_model("pocket_ref").atom})
    # resid_list.sort(key=lambda x: (x[2], int(x[1])))

    # print(f"\nОбщий карман, сформированный по всем HET-группам (радиус {radius} Å):")
    # for resn, resi, chain in resid_list:
    #     print(f"  - {resn:>3} {resi:>4} {chain}")

elif pocket_method == "2":
    failed_align = []  # список для белков, которые не выровнялись

    if init_align == "1":
        ref_sel = f"{ref_align_name} and chain {init_align_chain_id}"
        for name in names:
            if name == ref_align_name:
                continue
            mobil_str = f"{name} and chain {init_align_chain_id}"
            try:
                cmd.align(mobil_str, ref_sel)
                # aligned_path = os.path.join(f"{folder_path}/output/aligned_structures", f"{name}_{init_align_chain_id}_aligned_to_{ref_align_name}_{init_align_chain_id}.pdb")
                # cmd.save(aligned_path, mobil_str)
                # print(f"✅ {name}_{init_align_chain_id} сохранён как {name}_{init_align_chain_id}_aligned_to_{ref_align_name}_{init_align_chain_id}.pdb")
            except:
                print(f"❌ Не удалось выровнять {name}_{init_align_chain_id} по {ref_align_name}_{init_align_chain_id}")
                failed_align.append(name)  # Добавляем в список неудач

    elif init_align == "2":
        for name in names:
            if name == ref_align_name:
                continue
            mobil_str = f"{name}"
            try:
                cmd.align(mobil_str, ref_align_name)
                # aligned_path = os.path.join(f"{folder_path}/output/aligned_structures", f"{name}_aligned_to_{ref_align_name}.pdb")
                # cmd.save(aligned_path, mobil_str)
                # print(f"✅ {name} сохранён как {name}_aligned_to_{ref_align_name}.pdb")
            except:
                print(f"❌ Не удалось выровнять {name} по {ref_align_name}")
                failed_align.append(name)  # Добавляем в список неудач

    lig_sel = f"{ref_name} and resi {ligand_resi} and chain {ligand_chain}"
    cmd.select("ligand_ref", lig_sel)
    cmd.select("pocket_ref", f"(byres (ligand_ref around {radius})) and {ref_name} and polymer")
    resid_list = list({(atom.resn, atom.resi, atom.chain) for atom in cmd.get_model("pocket_ref").atom})
else:
    raise ValueError("Неверный метод задания кармана.")

resid_list.sort(key=lambda x: (int(x[1]), x[2]))


if pocket_method in ("1", "4", "5"):
    # --- Создание селекций карманов и Cα ---
    for name in names:
        parts = [f"(resi {resi} and chain {chain})" for resn, resi, chain in resid_list]
        pocket_sel = f"model {name} and polymer and (" + " or ".join(parts) + ")"
        ca_sel = f"{pocket_sel} and name CA"
        cmd.select(f"pocket_{name}", pocket_sel)
        cmd.select(f"pocket_CA_{name}", ca_sel)

elif pocket_method == "3":
    print("OK")

elif pocket_method == "2":
    for name in names:
        pocket_sel = f"(byres (ligand_ref around {radius})) and {name} and polymer"
        ca_sel = f"{pocket_sel} and name CA"
        cmd.select(f"pocket_{name}", pocket_sel)
        cmd.select(f"pocket_CA_{name}", ca_sel)

# --- Выравнивание всех по референсу ---
# ref_ca = f"pocket_CA_{ref_align_name}"
# os.makedirs(f"{folder_path}/output", exist_ok=True)
# os.makedirs(f"{folder_path}/output/aligned_structures", exist_ok=True)

# failed_align = []  # список для белков, которые не выровнялись
# for name in names:
#     if name == ref_align_name:
#         continue
#     moving_ca = f"pocket_CA_{name}"
#     try:
#         cmd.align(moving_ca, ref_ca)
#         aligned_path = os.path.join(f"{folder_path}/output/aligned_structures", f"{name}_aligned_to_{ref_align_name}.pdb")
#         cmd.save(aligned_path, name)
#         print(f"✅ {name} сохранён как {name}_aligned_to_{ref_align_name}.pdb")
#     except:
#         print(f"❌ Не удалось выровнять {name} по {ref_align_name}")
#         failed_align.append(name)  # Добавляем в список неудач



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

            if sel2_all == f"pocket_{ref_name}":
                aligned_path = os.path.join(f"{folder_path}/output/aligned_structures", f"{names[i]}_aligned_to_{ref_name}.pdb")
                cmd.save(aligned_path, sel1_all)
                print(f"✅ {names[i]} сохранён как {names[i]}_aligned_to_{ref_name}.pdb")


            rms_all_s = f"{rms_all:.5f}" if rms_all is not None else ""
            rms_ca_s = f"{rms_ca:.5f}" if rms_ca is not None else ""

            # print(f"{names[i]} vs {names[j]} | RMSD (all): {rms_all_s} | RMSD (Cα): {rms_ca_s}")
            
            rmsd_all.append((names[i], names[j], rms_all_s, reason_all))
            rmsd_ca.append((names[i], names[j], rms_ca_s, reason_ca))

            if rms_all is not None:
                heatmap_matrix[i, j] = rms_all
                heatmap_matrix[j, i] = rms_all  # симметрично

    # Удаление строк и столбцов, полностью заполненных NaN (кроме диагонали)
    def is_row_invalid(row, idx):
        return all(np.isnan(v) for j, v in enumerate(row) if j != idx)

    invalid_indices = [i for i, row in enumerate(heatmap_matrix) if is_row_invalid(row, i)]

    if invalid_indices:
        print("⚠️ Удалены структуры (выровнять не удалось):")
        for i in invalid_indices:
            print(f"  - {heatmap_labels[i]}")

        heatmap_matrix = np.delete(heatmap_matrix, invalid_indices, axis=0)
        heatmap_matrix = np.delete(heatmap_matrix, invalid_indices, axis=1)
        heatmap_labels = [name for i, name in enumerate(heatmap_labels) if i not in invalid_indices]
        n = heatmap_matrix.shape[0]
#----------
    # 1) Собираем DataFrame из матрицы для удобства
    df = pd.DataFrame(heatmap_matrix, index=heatmap_labels, columns=heatmap_labels)

    # 2) Заполняем NaN максимальным значением (чтобы не мешало clustering)
    max_val = np.nanmax(df.values)
    df_filled = df.fillna(max_val)

    # 3) linkage по строкам и столбцам
    row_Z = linkage(df_filled.values, method=linkage_method, metric='euclidean')
    col_Z = linkage(df_filled.values.T, method=linkage_method, metric='euclidean')

    # 4) Получаем порядок листьев
    row_order = leaves_list(row_Z)
    col_order = leaves_list(col_Z)

    # 5) Переставляем DataFrame
    df_ord = df_filled.iloc[row_order, :].iloc[:, col_order]
    ordered_labels = df_ord.index.tolist()

    # --- Добавляем средние значения ---
    df_ord['Mean'] = df_ord.mean(axis=1)       # среднее по строкам
    mean_row = df_ord.mean(axis=0)             # среднее по столбцам
    df_ord.loc['Mean'] = mean_row

    # --- Построение упорядоченной тепловой карты ---
    plt.figure(figsize=(12,10))
    sns.heatmap(
        df_ord,
        cmap='coolwarm',
        annot=True, fmt=".2f",
        annot_kws={"size": 4},
        linewidths=0.5,
        cbar_kws={"label":"RMSD (Å)"}
    )
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.title(f"Clustered RMSD heatmap ({linkage_method})")
    plt.tight_layout()

    heatmap_path = os.path.join(folder_path, "output", "rmsd_heatmap_clustered.png")
    plt.savefig(heatmap_path, dpi=300)
    plt.close()
    print(f"Упорядоченная тепловая карта сохранена как: rmsd_heatmap_clustered.png")

    # # Построение неупорядоченной тепловой карты
    # plt.figure(figsize=(12, 10))

    # # Вычисление среднего RMSD по строкам (игнорируя диагональ и nan)
    # mean_rmsd = []
    # for i in range(n):
    #     values = [heatmap_matrix[i, j] for j in range(n) if i != j and not np.isnan(heatmap_matrix[i, j])]
    #     avg = np.mean(values) if values else np.nan
    #     mean_rmsd.append(avg)

    # # Добавление среднего как последнего столбца
    # heatmap_matrix_with_mean = np.hstack((heatmap_matrix, np.array(mean_rmsd).reshape(-1, 1)))
    # heatmap_labels_with_mean = heatmap_labels + ["Mean"]

    # sns.heatmap(heatmap_matrix_with_mean, xticklabels=heatmap_labels_with_mean, yticklabels=heatmap_labels,
    #             cmap='coolwarm', annot=True, fmt=".2f", linewidths=0.5, cbar_kws={"label": "RMSD (Å)"})

    # plt.title(f"RMSD (all atoms) Heatmap – Method: {rmsd_method}")
    # plt.xticks(rotation=45, ha='right')
    # plt.yticks(rotation=0)
    # plt.tight_layout()

    # heatmap_path = os.path.join(folder_path, "output", "rmsd_heatmap.png")
    # plt.savefig(heatmap_path, dpi=300)
    # print(f"\nТепловая карта RMSD сохранена как: rmsd_heatmap.png")

    # Подготовка к Построению дендрограммы
    df = pd.DataFrame(heatmap_matrix, index=heatmap_labels, columns=heatmap_labels)


    # Удалим строки и столбцы, где больше 90% значений — NaN
    threshold = 0.7  # 70%
    print(df)
    row_nan_fraction = df.isna().mean(axis=1)
    col_nan_fraction = df.isna().mean(axis=0)

    df_filtered = df.loc[row_nan_fraction < threshold, df.columns[col_nan_fraction < threshold]]

    print("Максимальное значение RMSD:", df_filtered.values.max())
    print("Минимальное значение RMSD:", df_filtered.values.min())


    valid_values = df_filtered.mask(df_filtered == 0).stack()  # исключим диагональ (0.0)
    max_rmsd = valid_values.max()

    df_filtered.replace([np.inf, -np.inf], np.nan, inplace=True)

    nan_locs = np.argwhere(np.isnan(df_filtered.values))
    for i, j in nan_locs:
        row_name = df_filtered.index[i]
        col_name = df_filtered.columns[j]
        print(f"🔍 NaN между структурами: {row_name} ↔ {col_name}")

    df_filtered = df_filtered.fillna(max_rmsd)

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
            # Порог на основе 70% процента максимального расстояния слиянияя
            optimal_threshold = 0.7 * max(Z[:, 2])

            clusters = fcluster(Z, t=optimal_threshold, criterion='distance')
            clustering_description = f"knee_threshold_{optimal_threshold:.2f}"
            print(f"Порог, выбранный на основе 70% процента максимального расстояния слиянияя {optimal_threshold:.2f}")

        elif cl_choice == '4':
            # Автоматический выбор - Порог, выбранный методом локтя
            distances = Z[:, 2]
            kneedle = KneeLocator(range(1, len(distances)+1), distances, curve="convex", direction="increasing")
            optimal_threshold = distances[kneedle.knee]

            clusters = fcluster(Z, t=optimal_threshold, criterion='distance')
            clustering_description = f"knee_threshold_{optimal_threshold:.2f}"
            print(f"Порог, выбранный методом локтя (Elbow Method): {optimal_threshold:.2f}")
        else:
            print("Некорректный выбор.")

        plt.figure(figsize=(10, 7))
        dendrogram(Z, labels=df_filtered.index.tolist())
        plt.title(f"Иерархическая кластеризация на основе RMSD ({linkage_method})")
        plt.tight_layout()
        dendro_path = os.path.join(folder_path, "output", f"RMSD_dendrogram_{linkage_method}.png")
        plt.savefig(dendro_path)
        plt.close()
        print(f"Дендрограмма сохранена в: RMSD_dendrogram_{linkage_method}.png")

        if clusters is not None:
            cluster_df = pd.DataFrame({
                'Protein': df_filtered.index,
                'Cluster': clusters
            })
            output_csv_path = os.path.join(folder_path, "output", f"cluster_assignments_{clustering_description}_{linkage_method}.csv")
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
with open(f"{folder_path}/output/rmsd_all_atoms.csv", "w", newline='', encoding="utf-8-sig") as f:
    writer = csv.writer(f)
    writer.writerow(["Protein1", "Protein2", "RMSD_all_atoms", "Comment"])
    writer.writerows(rmsd_all)

with open(f"{folder_path}/output/rmsd_calpha.csv", "w", newline='', encoding="utf-8-sig") as f:
    writer = csv.writer(f)
    writer.writerow(["Protein1", "Protein2", "RMSD_Calpha", "Comment"])
    writer.writerows(rmsd_ca)

# --- Сохранение информации о карманах ---
with open(f"{folder_path}/output/info.txt", "w", encoding="utf-8-sig") as info_file:
    info_file.write(f"Структура для выравнивания: {ref_name}\n")
    if pocket_method == "1":
        info_file.write(f"Метод определения участка для выравнивания: На реф. структуре (радиус {radius} Å от остатка {ligand_resi} {ligand_chain})\n")
    elif pocket_method == "2":
        info_file.write(f"Метод определения участка для выравнивания: После предварительного выравнивания, для каждой структуры вокруг (радиус {radius} Å от остатка {ligand_resi} {ligand_chain})\n")
    elif pocket_method == "3":
        info_file.write(f"Метод определения участка для выравнивания: Для каждой структуры вокруг её HET-групп (в заданном {radius} Å радиусе)\n")
    elif pocket_method == "4":
        info_file.write(f"Метод определения участка для выравнивания: По введённому пользователем списку остатков (введённый список остатков) - {resi_chain} \n")
    elif pocket_method == "5":
        info_file.write(f"Метод определения участка для выравнивания: На реф. структуре по указанному идентификатору цепи {chain_id}\n")

    if pocket_method == '2':
        info_file.write(f"Референснуая структура для предварительного выравнивания: {ref_align_name}\n")
        if init_align == "1":
            info_file.write(f"Метод предварительного выравнивания: по цепи {init_align_chain_id}\n")
        elif init_align == "2":
            info_file.write(f"Метод предварительного выравнивания: по всей молекуле\n")

    info_file.write(f"метод расчёта RMSD: {rmsd_method}")

    info_file.write(f"\nВыбран метод кластеризации: {linkage_method}")

    if cl_choice == '1':
        info_file.write(f"\nВыбран Способ кластеризации: 1 — Задать количество кластеров (maxclust)")
        info_file.write(f"\nnum_clusters: {num_clusters}")
    elif cl_choice == '2':
        info_file.write(f"\nВыбран Способ кластеризации: 2 — Задать порог расстояния (distance)")
        info_file.write(f"\ndistance_threshold: {distance_threshold}")
    elif cl_choice == "3":
        info_file.write(f"\nВыбран Способ кластеризации: 3 — Автоматический: Порог на основе 70% процента максимального расстояния слиянияя") 
    elif cl_choice == "4":
        info_file.write(f"\nВыбран Способ кластеризации: 4 — Автоматический: Метод Изгиба или Локтя (Elbow Method)") 

    if comparison_mode == "1":
        info_file.write("\nВыбран Режим сравнения: 1 - Все со всеми (all vs all)")
    elif comparison_mode == "2":
        info_file.write("\nВыбран Режим сравнения: 2 - Все с референсом (all vs ref)")


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
        for atom in cmd.get_model(f"pocket_{ref_name}").atom
    }

    results = []

    # Собираем все diffs в список
    for name in names:
        if name == ref_name:
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

        results.append((name, diffs))

    # Сортируем: сначала те, у кого len(diffs)==0, затем по len(diffs) возрастанию
    sorted_results = sorted(
        results,
        key=lambda x: (len(x[1]) != 0, len(x[1]))
    )

# Пишем в файл уже в нужном порядке
    for name, diffs in sorted_results:
        if not diffs:
            info_file.write(f"{name}: идентичен референсу по карману.\n\n")
        else:
            info_file.write(f"{name}: отличается {len(diffs)} остатками:\n")
            for d in diffs:
                info_file.write(f"  {d}\n")
            info_file.write("\n")

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
        plt.savefig(f"{folder_path}/output/{filename}", dpi=300)
        plt.close()
        print(f"Гистограмма сохранена: {filename}")

    plot_rmsd_histogram(rmsd_all, "Гистограмма RMSD (все атомы)", "rmsd_all_atoms_hist.png")
    plot_rmsd_histogram(rmsd_ca, "Гистограмма RMSD (Cα)", "rmsd_calpha_hist.png")


# Сбор конфигурации
config = {
    "mode": mode,
    "folder_path": folder_path,
    "uniprot_id": uniprot_id,
    "pdb_id": pdb_id,
    "save_dir": save_dir,
    "do_preprocess": do_preprocess,
    "clean_options": to_remove,
    "ref_pocket": ref_name,
    "ref_align": ref_align_name,
    "init_align": init_align,
    "init_align_chain_id": init_align_chain_id,
    "pocket_method": pocket_method,
    "ligand_resi": ligand_resi,
    "ligand_chain": ligand_chain,
    "radius": radius,
    "resi_chain": resi_chain,
    "chain_id": chain_id,
    "comparison_mode": comparison_mode,
    "rmsd_method": rmsd_method_input,
    "linkage_method": linkage_method_choice,
    "cl_choice": cl_choice,
    "num_clusters": num_clusters,
    "distance_threshold": distance_threshold,
}



# Сохраняем
if not args.config:
    save_config_to_yaml(config, f"{folder_path}/output/run_config.yaml")


print("\n✅ Готово. Файлы сохранены:")
print("-----------------------------")
print(f"{folder_path}/output/")
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