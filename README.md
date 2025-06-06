# 🤖 PocketMaster 1.4 - Pocket Analyzer
PocketMaster — мощный и удобный инструмент для анализа сходства белковых карманов.
Скрипт помогает провести структурное выравнивание, анализировать карманы, вычислять RMSD, строить красивые графики и визуализации.

Данный скрипт предоставляет следующие возможности:

✅ Автоматически выполнять выравнивание структур белков

✅ Определять карман как на основе координат лиганда, так и вручную

✅ Рассчитывать RMSD по всем атомам и по атомам Cα в пределах карманов

✅ Поддержка нескольких методов RMSD (align, super, rms, rms_cur)

✅ Создавать информативные и наглядные графики и визуализации

✅ Сохранять выровненные структуры и результаты в формате CSV

✅ Генерировать текстовые отчёты и гистограммы с итогами анализа

# ⚙️ Зависимости

Для работы скрипта требуются следующие Python-библиотеки и программы:

-  Python 3.x
-  PyMOL (с поддержкой Python API, модуль pymol)
-  NumPy (`pip install numpy`)
-  Matplotlib (`pip install matplotlib`) — для построения гистограмм
-  Seaborn (`pip install seaborn`) — для улучшенной визуализации гистограмм


## 🚀 Как Запустить

Интерактивный PyMOL-скрипт помогает выполнять выравнивание белковых структур, строить карманы и рассчитывать RMSD.

Запустите скрипт с помощью Python:

	python PocketMaster1.4.py

  
## 📋 Как использовать

🧭 Следуйте интерактивным подсказкам:

🔹 **Укажите путь к папке с PDB-файлами**  
  (можно оставить пустым для текущей директории)

🔹 **Выберите референсные структуры:**  
  – для выравнивания  
  – для построения кармана

🔹 **Выберите метод задания кармана:**  
  – Радиальный (ввод остаток ID и радиус от него (Å))  
  – Ручной (ввод остатков вручную)  
  – По цепи (ввод chain ID)
 
🔹 **Укажите метод расчёта RMSD:**  
  – `align`строгая RMSD-оценка, используется для точного выравнивания, исключая атомы, не имеющие соответствий.  
  – `super` гибкий вариант, допускает несовпадения и автоматически подбирает соответствующие атомы.  
  – `rms` rms: точный RMSD, работает только при полном соответствии атомов (без несовпадений).  
  – `rms_cur` rms_cur: упрощённый и более быстрый расчёт RMSD, также требует полного соответствия атомов.  

🔹 **Выберите режим сравнения:**  
  – «все со всеми» (all vs all)  
  – «все с референсом» (all vs ref)  

🔹 **Выберите тип RMSD:**  
  – по всем атомам  
  – только по Cα-атомам

## После выполнения

Скрипт создаёт папку `aligned_output` в указанной директории, содержащую:

1. `aligned_structures/`: Выровненные PDB файлы (например, `structure_aligned_to_ref.pdb`).
2. `rmsd_all_atoms.csv`: RMSD для всех атомов кармана.
3. `rmsd_calpha.csv`: RMSD только для Cα атомов.
4. `info.txt`: Информация о референсной структуре, методе задания кармана,  выравниваниях и сравнении остатков карманов.
5. `rmsd_all_atoms_hist.png`: Гистограмма RMSD для всех атомов.
6. `rmsd_calpha_hist.png`: Гистограмма RMSD для Cα атомов.

## Структура директорий

- `data/` — папка с исходными PDB-файлами для анализа  
- `output/` — папка, куда после запуска скрипта автоматически сохраняются:  
  - выровненные структуры  
  - CSV-файлы с результатами расчетов  
  - текстовые отчёты и гистограммы

При запуске скрипта он автоматически обрабатывает все PDB-файлы из папки data/, выполняет выравнивание и анализ, после чего сохраняет все результаты в соответствующих форматах в папке data/output/. Такая организация упрощает работу с файлами и позволяет быстро ориентироваться в результатах


## 📌 Советы
Используй одинаково обработанные PDB-файлы (удаление воды, лигандов и т.д.)

## ⚠️ Возможные проблемы

- **Пустые селекции**  
  Если после определения кармана не удаётся выбрать атомы (например, из-за отсутствия соответствующих остатков), RMSD **не рассчитывается**.

- **Разное число атомов**  
  Скрипт выдаёт предупреждение при различии числа атомов в выравниваемых объектах. Это **особенно критично** для методов `rms` и `rms_cur`, которые требуют **полное соответствие атомов**.

- **Ошибки выравнивания**  
  Если структура не может быть выровнена (например, из-за отсутствующих атомов или несовпадения цепей), она **пропускается**, а информация об этом добавляется в **лог-файл**.

## 📧 Обратная связь
Если у тебя есть предложения или баг-репорты — не стесняйся связаться!

## ⭐ Поддержи проект
Если скрипт помог — поставь ⭐ на GitHub и расскажи другу-структурщику!

## Лицензия
MIT License. Используйте, модифицируйте и распространяйте свободно с указанием авторства.

