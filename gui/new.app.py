import streamlit as st
import os
import subprocess
import json
from pathlib import Path
from PIL import Image

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

# --- Tabs ---
tabs = st.tabs(["⚙️ Параметры", "📁 Результаты", "📦 Конфигурация"])

# ---------------- ПЕРВАЯ ВКЛАДКА ----------------
with tabs[0]:
    st.header("⚙️ Основные параметры")

    mode = st.selectbox("Выбери режим", ["1 — Локальная папка", "2 — UniProt ID", "3 — PDB ID"])
    folder_path = st.text_input("Путь к папке (если режим 1)", value=os.getcwd())
    uniprot_id = st.text_input("UniProt ID (если режим 2)")
    pdb_id = st.text_input("PDB ID (если режим 3)")

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

    aligned_path = Path(folder_path) / "aligned_output"

    if aligned_path.exists():
        show_file_if_exists(aligned_path / "rmsd_heatmap.png", "🧊 Тепловая карта RMSD")
        show_file_if_exists(aligned_path / f"RMSD_dendrogram_{config.get('linkage_method', 'ward')}.png", "🌳 Дендрограмма")
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
    uploaded_config = st.file_uploader("Загрузить JSON конфиг")

    if uploaded_config:
        config_data = json.load(uploaded_config)
        st.json(config_data)
