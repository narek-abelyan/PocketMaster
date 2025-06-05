import os
import csv
import numpy as np
import pandas as pd
from pymol import cmd
from scipy.cluster.hierarchy import dendrogram, linkage


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
        metals = ['ZN', 'MG', 'FE', 'CA', 'MN', 'CU', 'NA', 'K', 'CL']
        metal_selection = " or ".join([f"{name} and resn {ion}" for ion in metals])
        cmd.remove(f"{name} and solvent")
        cmd.remove(metal_selection)
    elif mode == "3":
        cmd.remove(f"{name} and not polymer")
    elif mode == "4":
        pass
    else:
        print("⚠️ Неверный режим очистки. Пропуск.")

# Использование
print_directory_contents_pretty(os.getcwd())

# Очистка среды
cmd.reinitialize()

# --- Ввод пути к папке с PDB ---
folder_path = input("\nВведи путь к папке с PDB файлами (по умолчанию текущая): ").strip()
if not folder_path:
    folder_path = os.getcwd()
files = [f for f in os.listdir(folder_path) if f.endswith('.pdb')]
files.sort()

print("\nДоступные PDB файлы:")
for i, f in enumerate(files, 1):
    print(f"{i}. {f}")

# --- Список ионов металлов ---
ions = ["MG", "ZN", "CA", "MN", "FE", "NA", "K", "CU", "CO", "NI"]

# --- Запрос на необходимость предварительной очистки ---
do_preprocess = input("\nПровести предварительную очистку структур? (1 - да, 0 - нет): ").strip()

to_remove = set()

if do_preprocess == '1':
    print("\nПожалуйста, выберите подходящие опции предварительной обработки структур (можно указать несколько пунктов последовательно):")
    print("1 — Удалить воду (solvent)")
    print("2 — Удалить ионы металлов (metal ions)")
    print("3 — Удалить всё, кроме белка (оставить только полимерную цепь)")
    print("4 — Удалить альтернативные конформации (altloc)")
    print("5 — Удалить атомы водорода (H)")
    print("6 — Удалить анизотропные параметры (ANISOU)")
    print("7 — Не очищать / Завершить выбор")
    print("8 - Сохранить обработанные структуры в определённой папке")
    

    save_dir = None  # Папка для сохранения обработанных структур
    to_remove = set()

    while True:
        choice = input("\nВведите номер опции (1–7), чтобы выбрать действие. или 7 — для завершения: ").strip()
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

print("\nДоступные PDB файлы:")
for i, f in enumerate(files, 1):
    print(f"{i}. {f}")

# --- Выбор двух референсных структур ---
print("\nВыберите референсную структуру для определения кармана:")
ref_index = int(input("Номер PDB-файла (карман): ")) - 1


print("\nВыберите референсную структуру для выравнивания:")
ref_align_index = int(input("Номер PDB-файла (выравнивание): ")) - 1

# Показать список всех HET-групп в выбранной структуре
ref_name = names[ref_index]
het_atoms = cmd.get_model(f"{ref_name} and not polymer and not solvent").atom
het_residues = sorted({(a.resn, a.resi, a.chain) for a in het_atoms}, key=lambda x: (int(x[1]), x[2]))

if het_residues:
    print(f" \nОбнаружены HET-группы в структуре {ref_name}:")
    for resn, resi, chain in het_residues:
        print(f"  - {resn:>3} {resi:>4} {chain}")
else:
    print(f"\nВ структуре {ref_name} не обнаружено HET-групп.")

# Выбор метода задания кармана
print(f"\nВыбери способ задания кармана:")
print("1 - Радиальное ( ввести остаток ID и радиус в Å)")
print("2 - Вручную (ввести список остатков)")
print("3 - По всей цепи (указать chain ID)")
print("4 - По всем HET-группам во всех структурах (в заданном радиусе Å)")
pocket_method = input("Введи способ (1, 2, 3 или 4): ").strip()

# Радиус окружности (по умолчанию 7)
ligand_resi = ""
radius = 0
resi_chain_input = ""
chain_id = ""

if pocket_method == "1":
    while True:
        ligand_input = input("\nВведи номер остатка и chain ID через пробел, например '40 A': ").strip()
        parts = ligand_input.split()
        if len(parts) != 2:
            print("Ошибка: нужно ввести номер остатка и chain ID через пробел, например '40 A'. Попробуй ещё раз.")
            continue
        ligand_resi, ligand_chain = parts
        ligand_chain = ligand_chain.upper()
        # Проверка, что ligand_resi — число
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

# Выбор режима сравнения
print("\nВыбери режим сравнения:")
print("1 - Все со всеми (all vs all)")
print("2 - Все с референсом (all vs ref) ")
comparison_mode = input("Введи режим (1 или 2): ").strip()

print("\nВыбери метод расчёта RMSD:")
print("1 - align   (Строгая оценка, исключает несоответствия)")
print("2 - super   (Гибкий подход, учитывает несовпадения)")
print("3 - rms     (Точное RMSD, требует полное соответствие атомов)")
print("4 - rms_cur (Быстрее, приближённое RMSD, тоже требует полное соответствие)")
rmsd_method_input = input("Введи метод (1, 2, 3 или 4): ").strip()

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


print("\n" + "-" * 30)
print("Запуск скрипта !")
print("-" * 30 + "\n")


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
color red, pocket_CA

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
ref_align_name = names[ref_align_index]
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
        print(f"✅ {name} сохранён как {ref_align_name}")
    except:
        print(f"❌ Не удалось выровнять {name} по {ref_align_name}")
        failed_align.append(name)  # Добавляем в список неудач


# --- Расчёт RMSD ---
rmsd_all = []
rmsd_ca = []

import matplotlib.pyplot as plt
import seaborn as sns

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

    # Убедиться, что размерность квадратная
    if df_filtered.shape[0] < 2:
        print("Недостаточно данных для кластеризации (менее двух структур).")
    else:
        # Преобразовать квадратную матрицу в конденсированную (condensed) форму
        from scipy.spatial.distance import squareform
        condensed = squareform(df_filtered.values, checks=False)

        # Построение дендрограммы
        Z = linkage(condensed, method='ward')
        plt.figure(figsize=(10, 7))
        dendrogram(Z, labels=df_filtered.index.tolist())
        plt.title("Иерархическая кластеризация на основе RMSD")
        plt.tight_layout()
        plt.savefig(os.path.join(folder_path, "aligned_output", "RMSD_dendrogram.png"))
        plt.close()
        print("Дендрограмма RMSD успешно построена: RMSD_dendrogram.png ")
        
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
print(f"- Гистограмма RMSD (Cα): rmsd_calpha_hist.png\n")
print("Дендрограмма RMSD успешно построена: RMSD_dendrogram.png ")
