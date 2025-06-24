import streamlit as st
import os
import subprocess
import json
from pathlib import Path
from PIL import Image

import subprocess
import sys

# try:
#     import py3Dmol
# except ImportError:
#     subprocess.check_call([sys.executable, "-m", "pip", "install", "py3dmol"])
#     import py3Dmol


# Импорт py3dmol для визуализации структур
try:
    import py3Dmol
    HAS_PY3DMOL = True
except ImportError:
    HAS_PY3DMOL = False

st.set_page_config(page_title="PocketMaster Web GUI", layout="wide")
st.title("🧬 PocketMaster 1.51 — Веб-интерфейс анализа карманов")

# --- Helper ---
def show_file_if_exists(path, caption=None, image=True):
    if os.path.exists(path):
        st.markdown(f"### {caption or os.path.basename(path)}")
        if image and path.endswith(('.png', '.jpg', '.jpeg')):
            st.image(Image.open(path), use_column_width=True)
        elif path.endswith(".csv"):
            import pandas as pd
            df = pd.read_csv(path)
            st.dataframe(df)
            st.download_button("⬇️ Скачать CSV", df.to_csv(index=False), file_name=os.path.basename(path))
        else:
            with open(path, 'r', encoding='utf-8') as f:
                content = f.read()
            st.code(content)

# Функция визуализации структуры из файла (PDB, MOL2, SDF)
def visualize_structure(filepath):
    if not HAS_PY3DMOL:
        st.warning("Для визуализации структур установите пакет py3dmol: pip install py3dmol")
        return

    with open(filepath, "r") as f:
        mol_str = f.read()

    ext = filepath.suffix.lower()
    viewer = py3Dmol.view(width=700, height=500)
    if ext == ".pdb":
        viewer.addModel(mol_str, "pdb")
    elif ext == ".mol2":
        viewer.addModel(mol_str, "mol2")
    elif ext == ".sdf":
        viewer.addModel(mol_str, "sdf")
    else:
        st.info(f"Визуализация для файлов {ext} не поддерживается")
        return
    viewer.setStyle({'stick': {}})
    viewer.zoomTo()
    viewer.show()

# --- Tabs ---
tabs = st.tabs(["⚙️ Параметры", "📁 Результаты", "📦 Конфигурация"])

# Глобальная переменная для конфигурации, чтобы переиспользовать в разных вкладках
if "config" not in st.session_state:
    st.session_state.config = {}

# ---------------- ПЕРВАЯ ВКЛАДКА ----------------
with tabs[0]:
    st.header("⚙️ Основные параметры")

    mode = st.selectbox("Выбери режим", ["1 — Локальная папка", "2 — UniProt ID", "3 — PDB ID"])
    folder_path = st.text_input("Путь к папке (если режим 1)", value=os.getcwd())
    uniprot_id = st.text_input("UniProt ID (если режим 2)")
    pdb_id = st.text_input("PDB ID (если режим 3)")

    # --- ВИЗУАЛИЗАЦИЯ СТРУКТУР ПРИ ВЫБОРЕ РЕЖИМА 1 ---
    if mode.startswith("1"):
        st.markdown("---")
        st.subheader("👓 Визуализация структур в папке")

        p = Path(folder_path)
        if p.exists() and p.is_dir():
            # Найдем подходящие файлы
            mol_files = list(p.glob("*.pdb")) + list(p.glob("*.mol2")) + list(p.glob("*.sdf"))
            if mol_files:
                selected_file = st.selectbox("Выберите структуру для визуализации", mol_files, format_func=lambda x: x.name)
                if selected_file:
                    visualize_structure(selected_file)
            else:
                st.info("В папке нет файлов структур с расширениями .pdb, .mol2 или .sdf")
        else:
            st.warning("Указанная папка не найдена или не является директорией.")

    st.markdown("---")
    st.subheader("🧽 Очистка")
    do_preprocess = st.checkbox("Включить очистку структур?", value=True)
    clean_options = st.multiselect("Опции очистки:", ["1", "2", "3", "4", "5", "6"])
    save_dir = st.text_input("Папка для сохранения (если требуется)", value="aligned_output")

    st.markdown("---")
    st.subheader("📌 Настройки кармана")
    pocket_method = st.selectbox("Метод задания кармана", ["1", "2", "3", "4"])
    ligand_resi = st.text_input("Номер остатка (для метода 1)", value="40")
    ligand_chain = st.text_input("Chain ID (для метода 1)", value="A")
    radius = st.number_input("Радиус (Å)", value=7.0)
    resi_chain = st.text_input("Список остатков (метод 2)", value="")
    chain_id = st.text_input("Chain ID (метод 3)", value="")

    st.markdown("---")
    st.subheader("📊 RMSD и кластеризация")
    comparison_mode = st.selectbox("Режим сравнения", ["1 — Все со всеми", "2 — Все с референсом"])
    rmsd_method = st.selectbox("Метод RMSD", ["1 — align", "2 — super", "3 — rms", "4 — rms_cur"])
    linkage_method = st.selectbox("Метод кластеризации", ["1 — ward", "2 — single", "3 — complete", "4 — average", "5 — centroid", "6 — median", "7 — weighted"])
    cl_choice = st.selectbox("Способ кластеризации", ["1 — по количеству", "2 — по порогу", "3 — авто"])
    num_clusters = st.number_input("Число кластеров", min_value=2, value=3)
    distance_threshold = st.number_input("Порог расстояния", value=3.0)

    st.markdown("---")
    run = st.button("🚀 Запустить анализ")

    if run:
        st.success("⏳ Запускаем анализ…")

        config = {
            "mode": int(mode[0]),
            "folder_path": folder_path,
            "uniprot_id": uniprot_id,
            "pdb_id": pdb_id,
            "do_preprocess": "1" if do_preprocess else "0",
            "clean_options": clean_options,
            "save_dir": save_dir,
            "pocket_method": pocket_method,
            "ligand_resi": ligand_resi,
            "ligand_chain": ligand_chain,
            "radius": radius,
            "resi_chain": resi_chain,
            "chain_id": chain_id,
            "comparison_mode": comparison_mode[0],
            "rmsd_method": rmsd_method[0],
            "linkage_method": linkage_method[0],
            "cl_choice": cl_choice[0],
            "num_clusters": int(num_clusters),
            "distance_threshold": float(distance_threshold)
        }

        st.session_state.config = config  # Сохраняем конфиг в сессию

        with open("temp_config.json", "w") as f:
            json.dump(config, f, indent=2)

        result = subprocess.run(
            ["python3", "PocketMaster1.51t.py", "--config", "temp_config.json"],
            capture_output=True, text=True
        )

        st.session_state.run_log = result.stdout + "\n" + result.stderr
        st.success("Готово! Перейдите во вкладку '📁 Результаты'")

# ---------------- ВТОРАЯ ВКЛАДКА ----------------
with tabs[1]:
    st.header("📁 Результаты анализа")

    if "run_log" in st.session_state:
        st.subheader("📋 Лог выполнения")
        st.code(st.session_state.run_log, language="bash")

    # Используем путь из конфигурации, если есть
    folder_for_results = st.session_state.config.get("folder_path", folder_path)
    save_dir = st.session_state.config.get("save_dir", "aligned_output")
    aligned_path = Path(folder_for_results) / save_dir

    if aligned_path.exists():
        show_file_if_exists(aligned_path / "rmsd_heatmap.png", "🧊 Тепловая карта RMSD")
        show_file_if_exists(aligned_path / f"RMSD_dendrogram_{st.session_state.config.get('linkage_method', 'ward')}.png", "🌳 Дендрограмма")
        show_file_if_exists(aligned_path / "rmsd_all_atoms_hist.png", "📉 Гистограмма RMSD (все атомы)")
        show_file_if_exists(aligned_path / "rmsd_calpha_hist.png", "📉 Гистограмма RMSD (Cα)")
        show_file_if_exists(aligned_path / "rmsd_all_atoms.csv", "📊 RMSD (все атомы)", image=False)
        show_file_if_exists(aligned_path / "rmsd_calpha.csv", "📊 RMSD (Cα)", image=False)
        show_file_if_exists(aligned_path / "info.txt", "📝 Информация", image=False)
    else:
        st.info("Результаты ещё не созданы или путь неверный.")

# ---------------- ТРЕТЬЯ ВКЛАДКА ----------------
with tabs[2]:
    st.header("📦 Работа с конфигурацией")

    uploaded_config = st.file_uploader("Загрузить JSON конфиг", type=["json"])

    if uploaded_config:
        try:
            config_data = json.load(uploaded_config)
            st.session_state.config = config_data
            st.success("Конфигурация загружена.")
            st.json(config_data)
        except Exception as e:
            st.error(f"Ошибка при загрузке конфигурации: {e}")

    if st.session_state.config:
        st.markdown("### Редактировать текущую конфигурацию")
        edited_config_str = st.text_area("JSON конфигурация", json.dumps(st.session_state.config, indent=2), height=300)

        if st.button("Сохранить конфигурацию в файл"):
            try:
                edited_config = json.loads(edited_config_str)
                with open("temp_config.json", "w") as f:
                    json.dump(edited_config, f, indent=2)
                st.success("Конфигурация сохранена в temp_config.json")
                st.session_state.config = edited_config
            except Exception as e:
                st.error(f"Ошибка при сохранении: {e}")
