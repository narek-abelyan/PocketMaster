from pymol import cmd
import tkinter as tk
from tkinter import filedialog, messagebox
import os
import csv

class PocketRMSDGUI:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("Pocket RMSD Calculator")
        self.folder_path = tk.StringVar()
        self.pocket_method = tk.StringVar(value="1")
        self.ligand_resi = tk.StringVar()
        self.radius = tk.DoubleVar(value=5.0)
        self.manual_residues = tk.StringVar()
        self.chain_id = tk.StringVar()
        self.comparison_mode = tk.StringVar(value="2")
        self.ref_index = tk.IntVar(value=1)
        self.files = []

        self.build_gui()

    def build_gui(self):
        tk.Label(self.root, text="Выбор папки:").grid(row=0, column=0, sticky="w")
        tk.Button(self.root, text="Выбрать", command=self.choose_folder).grid(row=0, column=1, sticky="w")

        self.files_listbox = tk.Listbox(self.root, height=5, width=50)
        self.files_listbox.grid(row=1, column=0, columnspan=2, pady=5)

        tk.Label(self.root, text="Метод определения кармана:").grid(row=2, column=0, sticky="w")
        tk.Radiobutton(self.root, text="По остатку лиганда", variable=self.pocket_method, value="1").grid(row=3, column=0, sticky="w")
        tk.Label(self.root, text="Номер остатка лиганда:").grid(row=4, column=0, sticky="w")
        tk.Entry(self.root, textvariable=self.ligand_resi).grid(row=4, column=1)

        tk.Label(self.root, text="Радиус вокруг лиганда:").grid(row=5, column=0, sticky="w")
        tk.Entry(self.root, textvariable=self.radius).grid(row=5, column=1)

        tk.Radiobutton(self.root, text="Ручной список остатков", variable=self.pocket_method, value="2").grid(row=6, column=0, sticky="w")
        tk.Label(self.root, text="Список остатков (через запятую):").grid(row=7, column=0, sticky="w")
        tk.Entry(self.root, textvariable=self.manual_residues).grid(row=7, column=1)

        tk.Radiobutton(self.root, text="Вся цепь", variable=self.pocket_method, value="3").grid(row=8, column=0, sticky="w")
        tk.Label(self.root, text="ID цепи:").grid(row=9, column=0, sticky="w")
        tk.Entry(self.root, textvariable=self.chain_id).grid(row=9, column=1)

        tk.Label(self.root, text="Режим сравнения:").grid(row=10, column=0, sticky="w")
        tk.Radiobutton(self.root, text="Все со всеми", variable=self.comparison_mode, value="1").grid(row=11, column=0, sticky="w")
        tk.Radiobutton(self.root, text="С одним референсом", variable=self.comparison_mode, value="2").grid(row=12, column=0, sticky="w")
        tk.Label(self.root, text="Индекс референса:").grid(row=13, column=0, sticky="w")
        tk.Entry(self.root, textvariable=self.ref_index).grid(row=13, column=1)

        tk.Button(self.root, text="Запустить анализ", command=self.run_analysis, bg="lightgreen").grid(row=14, column=0, columnspan=2, pady=10)

        self.root.mainloop()

    def choose_folder(self):
        folder = filedialog.askdirectory()
        if folder:
            self.folder_path.set(folder)
            self.files = [f for f in os.listdir(folder) if f.endswith(".pdb")]
            self.files_listbox.delete(0, tk.END)
            for f in self.files:
                self.files_listbox.insert(tk.END, f)

    def run_analysis(self):
        if not self.files:
            messagebox.showerror("Ошибка", "Сначала выберите папку с PDB-файлами.")
            return

        ref_idx = self.ref_index.get() - 1
        if not (0 <= ref_idx < len(self.files)):
            messagebox.showerror("Ошибка", "Неверный индекс референса.")
            return

        ref_name = os.path.splitext(self.files[ref_idx])[0]
        folder = self.folder_path.get()
        pocket_method = self.pocket_method.get()
        comp_mode = self.comparison_mode.get()

        cmd.reinitialize()
        objects = []

        for f in self.files:
            obj_name = os.path.splitext(f)[0]
            cmd.load(os.path.join(folder, f), obj_name)
            objects.append(obj_name)

        pockets = {}

        for obj in objects:
            if pocket_method == "1":
                sele = f"{obj} and resi {self.ligand_resi.get()}"
                pocket_sele = f"br. ({sele}) around {self.radius.get()}"
            elif pocket_method == "2":
                residues = self.manual_residues.get().replace(" ", "").split(",")
                pocket_sele = " or ".join([f"{obj} and resi {resi}" for resi in residues])
            elif pocket_method == "3":
                pocket_sele = f"{obj} and chain {self.chain_id.get()}"
            else:
                messagebox.showerror("Ошибка", "Неверный метод кармана.")
                return

            sele_name = f"pocket_{obj}"
            cmd.select(sele_name, pocket_sele)
            pockets[obj] = sele_name

        results = []
        comparisons = []

        if comp_mode == "2":
            ref = ref_name
            targets = [o for o in objects if o != ref]
            for target in targets:
                comparisons.append((ref, target))
        else:
            for i in range(len(objects)):
                for j in range(i + 1, len(objects)):
                    comparisons.append((objects[i], objects[j]))

        for ref, target in comparisons:
            rms_all = cmd.align(pockets[target], pockets[ref])[0]
            rms_ca = cmd.align(f"{pockets[target]} and name CA", f"{pockets[ref]} and name CA")[0]
            results.append((ref, target, round(rms_all, 3), round(rms_ca, 3)))

        out_file = os.path.join(folder, "rmsd_results.csv")
        with open(out_file, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Ref", "Target", "RMSD (all atoms)", "RMSD (Cα only)"])
            writer.writerows(results)

        messagebox.showinfo("Готово", f"Сравнение завершено.\nРезультаты сохранены в:\n{out_file}")

# Запуск GUI
PocketRMSDGUI()
