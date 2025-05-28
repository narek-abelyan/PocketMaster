import os
import csv
import numpy as np
from pymol import cmd
from tqdm import tqdm  # Импорт tqdm

def compute_rmsd(sel1, sel2):
    n1 = cmd.count_atoms(sel1)
    n2 = cmd.count_atoms(sel2)

    if n1 == 0 or n2 == 0:
        print(f"⚠️ Пропуск RMSD между {sel1} и {sel2} — пустая селекция.")
        return None, f"Пустая селекция: {n1} vs {n2}"
    
    allow_mismatch=int(n1*0.5)
    if abs(n1 - n2) > allow_mismatch:
        print(f"⚠️ Пропуск RMSD между {sel1} и {sel2} — слишком разное количество атомов ({n1} vs {n2}).")
        return None, f"Слишком Разное количество атомов: {n1} vs {n2}"
        
    try:
        rmsd = cmd.rms_cur(sel1, sel2)
        return rmsd, f"количество атомов ({n1} vs {n2})"
    except Exception as e:
        print(f"❌ Ошибка при расчёте RMSD между {sel1} и {sel2} --> количество атомов ({n1} vs {n2}): {e}")
        return None, str(e)

print(f"Текущая рабочая директория: {os.getcwd()}")

# Очистка среды
cmd.reinitialize()

# --- Ввод пользователя ---
folder_path = input("Введи путь к папке с PDB файлами (по умолчанию текущая): ").strip()
if not folder_path:
    folder_path = os.getcwd()
files = [f for f in os.listdir(folder_path) if f.endswith('.pdb')]
files.sort()

print("Доступные PDB файлы:")
for i, f in enumerate(files, 1):
    print(f"{i}. {f}")

ref_index = int(input("Введи номер референсного файла: ")) - 1

# --- Загрузка всех структур ---
names = []
for f in files:
    name = os.path.splitext(f)[0]
    names.append(name)
    cmd.load(os.path.join(folder_path, f), name)

# Показать список всех HET-групп в выбранной структуре
ref_name = names[ref_index]
het_atoms = cmd.get_model(f"{ref_name} and not polymer and not solvent").atom
het_residues = sorted({(a.resn, a.resi, a.chain) for a in het_atoms}, key=lambda x: (int(x[1]), x[2]))

if het_residues:
    print(f" Обнаружены HET-группы в структуре {ref_name}:")
    for resn, resi, chain in het_residues:
        print(f"  - {resn:>3} {resi:>4} {chain}")
else:
    print(f" В структуре {ref_name} не обнаружено HET-групп.")

# Выбор метода задания кармана
print(f"\nВыбери способ задания кармана:")
print("1 - Радиальное ( ввести остаток ID и радиус в Å)")
print("2 - Вручную (ввести список остатков)")
print("3 - По всей цепи (указать chain ID)")
pocket_method = input("Введи способ (1, 2 или 3): ").strip()

# Радиус окружности (по умолчанию 7)
ligand_resi = ""
radius = 0
resi_chain_input = ""
chain_id = ""

if pocket_method == "1":
    ligand_resi = input("Введи номер остатка лиганда в референсной структуре: ").strip()
    try:
        radius = float(input("Введи радиус для кармана (в Å, по умолчанию 7): ").strip())
    except ValueError:
        radius = 7.0
elif pocket_method == "2":
    resi_chain_input = input("Введи список остатков в формате (resn resi chain), например: PRO 46 A, ASN 61 A: ").strip()
elif pocket_method == "3":
    chain_id = input("Введи идентификатор цепи (chain ID), например: A: ").strip().upper()

# Выбор режима сравнения
print("Выбери режим сравнения:")
print("1 - Все со всеми (all vs all)")
print("2 - Все с референсом (all vs ref) ")
comparison_mode = input("Введи режим (1 или 2): ").strip()

print("\n" + "-" * 30)
print("Запуск скрипта !")
print("-" * 30 + "\n")


# --- Получение кармана ---
ref_name = names[ref_index]

if pocket_method == "1":
    lig_sel = f"{ref_name} and resi {ligand_resi}"
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

# --- Создание селекций карманов и Cα ---
for name in names:
    parts = [f"(resi {resi} and chain {chain})" for resn, resi, chain in resid_list]
    pocket_sel = f"{name} and polymer and (" + " or ".join(parts) + ")"
    ca_sel = f"{pocket_sel} and name CA"
    cmd.select(f"pocket_{name}", pocket_sel)
    cmd.select(f"pocketCA_{name}", ca_sel)

# --- Выравнивание всех по референсу ---
ref_ca = f"pocketCA_{ref_name}"
os.makedirs("aligned_structures", exist_ok=True)

failed_align = []  # список для белков, которые не выровнялись

print("Выравнивание структур:")
for name in tqdm(names, desc="Aligning", unit="structure"):
    if name == ref_name:
        continue
    moving_ca = f"pocketCA_{name}"
    try:
        cmd.align(moving_ca, ref_ca)
        aligned_path = os.path.join("aligned_structures", f"{name}_aligned_to_{ref_name}.pdb")
        cmd.save(aligned_path, name)
        print(f"✅ {name} сохранён как {aligned_path}")
    except:
        print(f"❌ Не удалось выровнять {name} по {ref_name}")
        failed_align.append(name)  # Добавляем в список неудач


# --- Расчёт RMSD ---
rmsd_all = []
rmsd_ca = []

print("Расчёт RMSD:")
if comparison_mode == "1":
    # All vs all
    total = len(names) * (len(names) - 1) // 2
    with tqdm(total=total, desc="RMSD all vs all", unit="pair") as pbar:
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                n1 = names[i]
                n2 = names[j]
                sel1_all = f"pocket_{n1}"
                sel2_all = f"pocket_{n2}"
                sel1_ca = f"pocketCA_{n1}"
                sel2_ca = f"pocketCA_{n2}"

                rms_all, reason_all = compute_rmsd(sel1_all, sel2_all)
                rms_ca_val, reason_ca = compute_rmsd(sel1_ca, sel2_ca)
                
                rms_all_s = f"{rms_all:.5f}" if rms_all is not None else ""
                rms_ca_s = f"{rms_ca_val:.5f}" if rms_ca_val is not None else ""

                rmsd_all.append((n1, n2, rms_all_s, reason_all))
                rmsd_ca.append((n1, n2, rms_ca_s, reason_ca))
                pbar.update(1)

elif comparison_mode == "2":
    # Все с референсом
    with tqdm(total=len(names) - 1, desc="RMSD all vs ref", unit="structure") as pbar:
        for name in names:
            if name == ref_name:
                continue
            sel1_all = f"pocket_{ref_name}"
            sel2_all = f"pocket_{name}"
            sel1_ca = f"pocketCA_{ref_name}"
            sel2_ca = f"pocketCA_{name}"

            rms_all, reason_all = compute_rmsd(sel1_all, sel2_all)
            rms_ca_val, reason_ca = compute_rmsd(sel1_ca, sel2_ca)
            
            rms_all_s = f"{rms_all:.5f}" if rms_all is not None else ""
            rms_ca_s = f"{rms_ca_val:.5f}" if rms_ca_val is not None else ""

            rmsd_all.append((ref_name, name, rms_all_s, reason_all))
            rmsd_ca.append((ref_name, name, rms_ca_s, reason_ca))
            pbar.update(1)
else:
    raise ValueError("Неверный режим сравнения.")


# --- Сохранение RMSD ---
with open("rmsd_all_atoms.csv", "w", newline='') as f_all:
    writer = csv.writer(f_all)
    writer.writerow(["Structure 1", "Structure 2", "RMSD (all atoms)", "Примечание"])
    writer.writerows(rmsd_all)

with open("rmsd_ca_atoms.csv", "w", newline='') as f_ca:
    writer = csv.writer(f_ca)
    writer.writerow(["Structure 1", "Structure 2", "RMSD (Cα atoms)", "Примечание"])
    writer.writerows(rmsd_ca)

print("\n✅ RMSD расчёты сохранены в файлы 'rmsd_all_atoms.csv' и 'rmsd_ca_atoms.csv'.")
