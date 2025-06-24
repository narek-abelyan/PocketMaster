import streamlit as st
import os
import subprocess

st.set_page_config(layout="wide")
st.title("🧪 PocketMaster 1.51 — Веб-интерфейс")

# --- Выбор режима ---
mode = st.selectbox("Выбери режим работы", ["Локальная папка (1)", "UniProt ID (2)", "PDB ID → UniProt (3)"])
mode_map = {"Локальная папка (1)": "1", "UniProt ID (2)": "2", "PDB ID → UniProt (3)": "3"}
mode_val = mode_map[mode]

# --- Ввод параметров ---
folder_path = st.text_input("Путь к папке с PDB (для режима 1)", value=os.getcwd())
uniprot_id = st.text_input("UniProt ID (для режима 2)", value="")
pdb_id = st.text_input("PDB ID (для режима 3)", value="")

config_path = st.file_uploader("Загрузить JSON-конфигурацию", type=["json"])

run_button = st.button("🚀 Запустить анализ")

# --- Запуск ---
if run_button:
    with st.spinner("Запуск скрипта..."):

        config_file = "temp_config.json"

        # если загружен конфиг — сохранить
        if config_path is not None:
            with open(config_file, "wb") as f:
                f.write(config_path.read())
        else:
            # иначе — сгенерировать минимальный config
            import json
            config_data = {
                "mode": int(mode_val),
                "folder_path": folder_path if mode_val == "1" else None,
                "uniprot_id": uniprot_id if mode_val == "2" else None,
                "pdb_id": pdb_id if mode_val == "3" else None
            }
            with open(config_file, "w") as f:
                json.dump(config_data, f, indent=2)

        # вызов оригинального скрипта
        result = subprocess.run(
            ["python", "PocketMaster1.51t.py", "--config", config_file],
            capture_output=True,
            text=True
        )

        st.subheader("Лог выполнения")
        st.code(result.stdout + "\n" + result.stderr)

        if result.returncode == 0:
            st.success("Анализ завершён ✅")
            st.markdown(f"📁 Проверь папку: `{folder_path}/aligned_output/`")
        else:
            st.error("❌ Ошибка при выполнении скрипта.")
