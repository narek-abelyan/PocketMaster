import streamlit as st
import os
import json
import subprocess
from pathlib import Path
from PIL import Image
import streamlit.components.v1 as components

# Проверка наличия py3Dmol
try:
    import py3Dmol
    HAS_PY3DMOL = True
except ImportError:
    HAS_PY3DMOL = False

st.set_page_config(page_title="PocketMaster Web GUI", layout="wide")
st.title("🧬 PocketMaster 1.51 — Веб-интерфейс анализа карманов")

# --- Helper ---
def show_file_if_exists(path, caption=None, image=True):
    path = Path(path)
    if path.exists():
        st.markdown(f"### {caption or path.name}")
        if image and path.suffix.lower() in (".png", ".jpg", ".jpeg"):
            st.image(Image.open(path), use_column_width=True)
        elif path.suffix.lower() == ".csv":
            import pandas as pd
            df = pd.read_csv(path)
            st.dataframe(df)
            st.download_button("⬇️ Скачать CSV", df.to_csv(index=False), file_name=path.name)
        else:
            with open(path, "r", encoding="utf-8") as f:
                content = f.read()
            st.code(content)

# --- Интерактивная визуализация структуры ---
def visualize_structure_interactive(filepath):
    if not HAS_PY3DMOL:
        st.warning("Установите py3Dmol: pip install py3dmol")
        return

    filepath = Path(filepath)
    if not filepath.exists():
        st.warning(f"Файл {filepath} не найден")
        return

    with open(filepath, "r") as f:
        mol_str = f.read()

    ext = filepath.suffix.lower()
    format_map = {".pdb": "pdb", ".mol2": "mol2", ".sdf": "sdf"}
    fmt = format_map.get(ext)
    if not fmt:
        st.info(f"Формат {ext} не поддерживается для визуализации")
        return

    # Интерфейс управления стилями
    col1, col2, col3 = st.columns(3)
    with col1:
        style = st.selectbox("Стиль", ["stick", "line", "cartoon", "surface", "sphere", "cross"], index=0)
    with col2:
        color = st.selectbox("Цвет", ["spectrum", "white", "red", "green", "blue", "yellow", "cyan", "magenta"], index=0)
    with col3:
        opacity = st.slider("Прозрачность (только для surface)", 0.0, 1.0, 0.6, 0.05)

    # Выбор цепей (если есть в файле)
    chain_input = st.text_input("Вывести только цепи (через запятую, например A,B,C)", value="")
    chains = [c.strip() for c in chain_input.split(",") if c.strip()]

    viewer = py3Dmol.view(width=850, height=600)
    viewer.addModel(mol_str, fmt)

    style_dict = {}
    if style == "surface":
        style_dict = {style: {"opacity": opacity, "color": color}}
    else:
        style_dict = {style: {"color": color}}

    if chains:
        for chain in chains:
            viewer.setStyle({"chain": chain}, style_dict)
    else:
        viewer.setStyle({}, style_dict)

    viewer.zoomTo()
    viewer.setBackgroundColor("0xeeeeee")
    viewer.rotate(90, "y")

    html = viewer._make_html()
    components.html(html, height=600, scrolling=False)

# --- Tabs ---
tabs = st.tabs(["⚙️ Параметры", "📁 Результаты", "📦 Конфигурация"])

if "config" not in st.session_state:
    st.session_state.config = {}

# ---------------- ПЕРВАЯ ВКЛАДКА ----------------
with tabs[0]:
    st.header("⚙️ Основные параметры")

    mode = st.selectbox("Выбери режим", ["1 — Локальная папка", "2 — UniProt ID", "3 — PDB ID"])
    folder_path = st.text_input("Путь к папке (если режим 1)", value=os.getcwd())
    uniprot_id = st.text_input("UniProt ID (если режим 2)")
    pdb_id = st.text_input("PDB ID (если режим 3)")

    if mode.startswith("1"):
        st.subheader("👓 Визуализация структур")
        p = Path(folder_path)
        if p.exists():
            files = list(p.glob("*.pdb")) + list(p.glob("*.mol2")) + list(p.glob("*.sdf"))
            if files:
                selected_file = st.selectbox("Файл структуры", files, format_func=lambda x: x.name)
                if selected_file:
                    visualize_structure_interactive(selected_file)
            else:
                st.info("В папке нет поддерживаемых файлов.")
        else:
            st.warning("Папка не найдена.")

    st.subheader("🧽 Очистка")
    do_preprocess = st.checkbox("Включить очистку структур?", value=True)
    clean_options = st.multiselect("Опции очистки", ["1", "2", "3", "4", "5", "6"])
    save_dir = st.text_input("Папка для сохранения", value="aligned_output")

    st.subheader("📌 Настройки кармана")
    pocket_method = st.selectbox("Метод задания кармана", ["1", "2", "3", "4"])
    ligand_resi = st.text_input("Номер остатка", value="40")
    ligand_chain = st.text_input("Chain ID", value="A")
    radius = st.number_input("Радиус", value=7.0)
    resi_chain = st.text_input("Список остатков", value="")
    chain_id = st.text_input("Chain ID (метод 3)", value="")

    st.subheader("📊 RMSD и кластеризация")
    comparison_mode = st.selectbox("Режим сравнения", ["1 — Все со всеми", "2 — Все с референсом"])
    rmsd_method = st.selectbox("Метод RMSD", ["1 — align", "2 — super", "3 — rms", "4 — rms_cur"])
    linkage_method = st.selectbox("Метод кластеризации", ["1 — ward", "2 — single", "3 — complete", "4 — average", "5 — centroid", "6 — median", "7 — weighted"])
    cl_choice = st.selectbox("Тип кластеризации", ["1 — по количеству", "2 — по порогу", "3 — авто"])
    num_clusters = st.number_input("Число кластеров", min_value=2, value=3)
    distance_threshold = st.number_input("Порог расстояния", value=3.0)

    if st.button("🚀 Запустить анализ"):
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
        st.session_state.config = config
        with open("temp_config.json", "w") as f:
            json.dump(config, f, indent=2)
        result = subprocess.run(["python3", "PocketMaster1.51t.py", "--config", "temp_config.json"], capture_output=True, text=True)
        st.session_state.run_log = result.stdout + "\n" + result.stderr
        st.success("Готово! Перейдите во вкладку '📁 Результаты'")

# ---------------- ВТОРАЯ ВКЛАДКА ----------------
with tabs[1]:
    st.header("📁 Результаты анализа")

    if "run_log" in st.session_state:
        st.subheader("📋 Лог")
        st.code(st.session_state.run_log, language="bash")

    config = st.session_state.get("config", {})
    result_dir = Path(config.get("folder_path", ".")) / config.get("save_dir", "aligned_output")

    if result_dir.exists():
        show_file_if_exists(result_dir / "rmsd_heatmap.png", "🧊 Тепловая карта RMSD")
        show_file_if_exists(result_dir / f"RMSD_dendrogram_{config.get('linkage_method', 'ward')}.png", "🌳 Дендрограмма")
        show_file_if_exists(result_dir / "rmsd_all_atoms_hist.png", "📉 RMSD (все атомы)")
        show_file_if_exists(result_dir / "rmsd_calpha_hist.png", "📉 RMSD (Cα)")
        show_file_if_exists(result_dir / "rmsd_all_atoms.csv", "📊 RMSD (все атомы)", image=False)
        show_file_if_exists(result_dir / "rmsd_calpha.csv", "📊 RMSD (Cα)", image=False)
        show_file_if_exists(result_dir / "info.txt", "📝 Информация", image=False)
    else:
        st.info("Результаты не найдены")

# ---------------- ТРЕТЬЯ ВКЛАДКА ----------------
with tabs[2]:
    st.header("📦 Работа с конфигурацией")

    uploaded = st.file_uploader("Загрузить конфиг JSON", type=["json"])
    if uploaded:
        try:
            config_data = json.load(uploaded)
            st.session_state.config = config_data
            st.success("Конфигурация загружена")
            st.json(config_data)
        except Exception as e:
            st.error(f"Ошибка загрузки: {e}")

    if st.session_state.get("config"):
        st.markdown("### Редактирование текущей конфигурации")
        text = st.text_area("Редактировать JSON", json.dumps(st.session_state.config, indent=2))
        if st.button("💾 Сохранить конфигурацию"):
            try:
                parsed = json.loads(text)
                with open("temp_config.json", "w") as f:
                    json.dump(parsed, f, indent=2)
                st.success("Конфигурация сохранена")
            except Exception as e:
                st.error(f"Ошибка сохранения: {e}")
