## EN

# 🤖 PocketMaster — Pocket Analyzer

PocketMaster is a flexible and automated tool for analyzing, clustering, and visualizing protein binding sites. It allows you to quickly compare structures, explore functional regions of proteins, and generate clear results, even when working with hundreds or thousands of models. This makes it particularly useful in the early stages of drug design, when analyzing and selecting the correct protein structures is crucial.

---

## Main Functionalities of PocketMaster

* ✅ **Automatic structure alignment**
* ✅ **Flexible methods for defining binding sites**
* ✅ **Support for multiple RMSD methods:** `align`, `cealign`, `super`, `rms`, `rms_cur`
* ✅ **RMSD calculation across all atoms or only Cα atoms within the binding sites**
* ✅ **Creation of clear visualizations: heatmaps, dendrograms, histograms**
* ✅ **Saving of aligned structures** and analysis results
* ✅ **Generation of summary reports with structural and sequence differences of binding sites
* ✅ **Choice of clustering linkage methods:** `ward`, `single`, `complete`, `average`, `centroid`, `median`, `weighted`
* ✅ **Customizable clustering options: by number of clusters (maxclust), by distance threshold (distance), or automatic modes (Elbow method, 70% of maximum linkage distance threshold)**

With PocketMaster, you can quickly analyze the similarity of protein binding sites and obtain detailed visualizations and quantitative assessments for in-depth structural insights.

---

## Visualization Examples

<p float="left">
  <img src="images/rmsd_heatmap.png" height="280" />
  <img src="images/rmsd_all_atoms_hist.png" height="220" />
  <img src="images/RMSD_dendrogram_average.png" height="220" />
</p>

**Pocket comparison visualizations:** `heatmap.png`, `rmsd_all_atoms_hist.png`, `dendrogram.png`

---

## ⚙️ Dependencies

### Required Libraries and Tools

The following libraries and tools are required for the script to run properly:

- **Python 3.x**
- **PyMOL** with Python API support (module `pymol`) — *install separately*
- **NumPy** — numerical computations  
  Installation: `pip install numpy`
- **Matplotlib** — plotting graphs and histograms  
  Installation: `pip install matplotlib`
- **Seaborn** — enhanced visualizations (including heatmaps)  
  Installation: `pip install seaborn`
- **Pandas** — working with tables and CSV files  
  Installation: `pip install pandas`
- **SciPy** — for clustering and distance analysis  
  Installation: `pip install scipy`
- **YAML** — working with configuration YAML files  
  Installation: `pip install yaml`

> ⚠️ Note: PyMOL is **not** installed via `pip`. Install PyMOL separately (manually or via `conda`).

---

## Installation

You can install the project dependencies in two convenient ways:

### 1. Install via `pip` (using `requirements.txt`)

If you already have Python and `pip`, run:

```bash
pip install -r requirements.txt
```

> ⚠️ Note: PyMOL cannot be installed via pip.  
> It must be installed separately — either manually or via conda (see below).

### 2. Install via `conda` (using `environment.yml`) - recommended

Recommended method: use Anaconda or Miniconda. Create and activate the environment with:

```bash
conda env create -f environment.yml
conda activate pmaster
```

---

## 🚀 How to run

Two running modes are supported:

### Interactive Mode

Run the script interactively, entering parameters step by step:

```bash
python PocketMaster.py
```

### Configuration File Mode

Run the script automatically using a YAML configuration file with pre-defined parameters:

```bash
python PocketMaster.py --config config.yaml
```

📌 See an example configuration in `examples/config.yaml`.

When running in interactive mode, a run_config.txt file is generated at the end. It automatically saves all parameters entered by the user, allowing you to run the script later in automatic mode with the same settings, without re-entering parameters.

---

## 📋 How to Use

Run the script **PocketMaster.py** in interactive mode using the command shown above  
to enter interactive mode and set the parameters manually.

🧭 Follow the interactive prompts:

🔹 **Select the operating mode:**  
  1 – Use a local folder with PDB files  
  2 – Download all corresponding PDB structures based on a UniProt ID  
  3 – Determine the UniProt ID from a PDB ID and download all corresponding PDB structures

🔹 **Perform preliminary structure cleanup?**  
  1 – Yes  
  2 – No  

🔹 **Select the appropriate structure preprocessing options:**  
  1 – Remove water (solvent)  
  2 – Remove ions (Cu, CL, etc.)  
  3 – Remove sulfates and phosphates (SO4, PO4, etc.)  
  4 – Remove buffer components (TRS, MES, HEP, etc.)  
  5 – Remove cryoprotectants (GOL, EDO, MPD, etc.)  
  6 – Remove reducing agents (DTT, BME, TCEP)  
  7 – Remove all water, ions, buffers, cryoprotectants, phosphates, and reducing agents  
  8 – Remove modified amino acid residues (CSO, MSE, SEP, TPO, PTR, etc.)  
  9 – Remove everything except the protein (keep only the polymer chain)  
  10 – Remove alternate conformations (altloc)  
  11 – Remove anisotropic parameters (ANISOU)  
  12 – Remove hydrogen atoms (H)  
  13 – Save processed structures to a specified folder  
  14 – Do not clean / Finish selection  

🔹 **Select the reference structure for alignment:**  
  – Enter its number from the provided list  

🔹 **Choose the method for defining the alignment region:**  
  1 – On the reference structure using a specified residue ID and radius (Å), then search and align in all structures  
  2 – After preliminary alignment of all structures with each other, define the region around the selected reference residue for each structure  
  3 – For each structure, around its HET groups within the specified radius (Å)  
  4 – Using a user-provided list of residues, then search in all structures  
  5 – On the reference structure using a specified chain ID, then search and align in all structures

---

<p align="left">Supported methods for defining alignment regions and their detailed interpretation:</p>
<p align="left">
  <img src="images/graphic.png" height="400" />
</p>

---

🔹 **Select the comparison mode:**  
  1 – "all vs all"  
  2 – "all vs reference"  

🔹 **Specify the RMSD calculation method:**  
  1 – `align`: strict RMSD evaluation, used for precise alignment, excluding unmatched atoms.  
  2 – `cealign`: geometric-based structural alignment, effective even with low similarity.  
  3 – `super`: flexible option, allows mismatches and automatically matches corresponding atoms.  
  4 – `rms`: accurate RMSD, works only with complete atom correspondence (no mismatches).  
  5 – `rms_cur`: simplified and faster RMSD calculation, also requires full atom correspondence.  

🔹 **Select the hierarchical clustering method:**  
  1 – ward (minimizes intra-cluster variance, requires Euclidean distance)  
  2 – single (minimum distance between clusters)  
  3 – complete (maximum distance between clusters)  
  4 – average (average distance between clusters, UPGMA)  
  5 – centroid (distance between cluster centroids)  
  6 – median (median distance between clusters)  
  7 – weighted (weighted average distance, WPGMA)


---

<h3 align="left">📊 Comparison of Clustering Methods</h3>

<p align="left">

| Method      | Compact Clusters   | Elongated Clusters  | Sensitivity to Outliers |
|------------|------------------|------------------|------------------------|
| Ward       | ✅ Excellent       | ❌ Poor           | ⚠️ Medium             |
| Single     | ❌ Poor           | ✅ Excellent      | ⚠️ High               |
| Complete   | ✅ Good           | ❌ Poor           | ✅ Stable             |
| Average    | ✅ Versatile      | ✅ Fairly Good    | ⚠️ Medium             |
| Centroid   | ⚠️ Unstable       | ⚠️ Unstable      | ⚠️ Unstable           |
| Median     | ⚠️ Unstable       | ⚠️ Unstable      | ⚠️ Unstable           |
| Weighted   | ✅ Good           | ✅ Good          | ⚠️ Medium             |

</p>
<br>

🔹 **Clustering parameters: number of clusters or distance threshold**  
  1 – Specify the number of clusters (`maxclust`)  
  2 – Specify a distance threshold (`distance`)  
  3 – Automatic: threshold based on 70% of the maximum merge distance  
  4 – Automatic: Elbow Method  

🔹 **Select the RMSD type:**  
  1 – All atoms  
  2 – Only Cα atoms

---

## After Execution

The script creates an `output` folder in the specified directory, containing:

1. `aligned_structures/`: Aligned PDB files (e.g., `structure_aligned_to_ref.pdb`).  
2. `rmsd_all_atoms.csv`: RMSD for all pocket atoms.  
3. `rmsd_calpha.csv`: RMSD for Cα atoms only.  
4. `info.txt`: Information about the reference structure, pocket definition method, alignments, and pocket residue comparisons.  
5. `rmsd_all_atoms_hist.png`: RMSD histogram for all atoms.  
6. `rmsd_calpha_hist.png`: RMSD histogram for Cα atoms.  
7. `Clustering results`: `cluster_assignments.csv`.  
8. `Heatmap`: `rmsd_heatmap.png`.  
9. `Dendrogram`: `RMSD_dendrogram.png`.  

## Directory Structure

- `data/` — folder containing the original PDB files for analysis.  
- `output/` — folder where the script automatically saves after execution:  
  - aligned structures  
  - CSV files with calculation results  
  - text reports and histograms  

When run, the script automatically processes all PDB files in the `data/` folder, performs alignment and analysis, and saves all results in the corresponding formats inside `data/output/`. This organization simplifies file management and allows you to quickly navigate the results.

---

## 📌 Tips
Use consistently preprocessed PDB files (e.g., with water, ligands, etc. removed).

## ⚠️ Possible Issues

- **Empty selections**  
  If no atoms can be selected after defining the pocket (e.g., due to missing residues), RMSD **will not be calculated**.

- **Different number of atoms**  
  The script issues a warning if the number of atoms in the structures being aligned differs. This is **especially critical** for the `rms` and `rms_cur` methods, which require **full atom correspondence**.

- **Alignment errors**  
  If a structure cannot be aligned (e.g., due to missing atoms or chain mismatches), it is **skipped**, and information about the issue is added to the **log file**.

## 📧 Feedback
If you have suggestions or bug reports, feel free to get in touch! (narek.abelyan@gmail.com)

## ⭐ Support the Project
If the script was helpful — give it a ⭐ on GitHub and share it with a fellow structural biologist!

## License
MIT License. Feel free to use, modify, and distribute with proper attribution.

***
## RU

# 🤖 PocketMaster - Pocket Analyzer
PocketMaster - удобный инструмент для анализа сходства белковых карманов.
Скрипт помогает провести структурное выравнивание, анализировать карманы, вычислять RMSD, строить красивые графики и визуализации.

## Основные возможности

- ✅ **Автоматическое выравнивание** белковых структур  
- ✅ **Гибкое определение кармана** - как по координатам лиганда, так и вручную  
- ✅ **Расчёт RMSD** по всем атомам и по атомам Cα в пределах карманов  
- ✅ **Поддержка различных методов RMSD**: `align`, `cealign`, `super`, `rms`, `rms_cur`  
- ✅ **Сохранение выровненных структур** и результатов анализа  
- ✅ **Генерация кратких отчётов** с данными о структурных и последовательностных различиях карманов  
- ✅ **Построение наглядных визуализаций**: Тепловые карты (heatmap), Дендрограммы,Гистограммы  
- ✅ **Выбор метода кластеризации**: `ward`, `single`, `complete`, `average`, `centroid`, `median`, `weighted`
- ✅ **Настройка способа кластеризации**: по числу кластеров (`maxclust`), по порогу расстояния (`distance`), автоматические режимы (`Elbow method`, `70% of maximum linkage distance threshold`)

---

С PocketMaster вы сможете быстро и эффективно анализировать пространственное сходство активных участков белков, получая подробные отчёты и визуализации для глубокой интерпретации результатов.

## Примеры визуализаций

<p float="left">
  <img src="images/rmsd_heatmap.png" height="280" />
  <img src="images/rmsd_all_atoms_hist.png" height="215" />
  <img src="images/RMSD_dendrogram_average.png" height="215" />
</p>

**Сравнение карманов:** heatmap.png, rmsd_all_atoms_hist.png, dendrogram.png

## ⚙️ Зависимости

Для корректной работы скрипта необходимы следующие библиотеки и инструменты:

- **Python 3.x**
- **PyMOL** с поддержкой Python API (доступ к модулю `pymol`)
- **NumPy** — численные вычисления  
  Установка: `pip install numpy`
- **Matplotlib** — построение графиков и гистограмм  
  Установка: `pip install matplotlib`
- **Seaborn** — улучшенная визуализация (в том числе тепловые карты)  
  Установка: `pip install seaborn`
- **Pandas** — работа с таблицами и CSV  
  Установка: `pip install pandas`
- **SciPy** — для кластеризации и анализа расстояний  
  Установка: `pip install scipy`
- **YAML** — для работы с конфигурационными YAML-файлами  
  Установка: `pip install yaml`

## Установка зависимостей

Для удобства вы можете установить все необходимые зависимости двумя способами:

### 1. Установка через pip (с использованием файла requirements.txt)

Если у вас уже установлен Python и менеджер пакетов pip, выполните команду:

```bash
pip install -r requirements.txt
```
> ⚠️ Обратите внимание: PyMOL не устанавливается через pip.  
> Его необходимо устанавливать отдельно — вручную или через conda (см. ниже).

### 2. Установка через conda (с использованием environment.yml)

Рекомендуемый способ - использовать Anaconda или Miniconda. Для создания и активации окружения выполните:

```bash
conda env create -f environment.yml
conda activate pmaster
```

## 🚀 Как запустить

Поддерживаются два режима запуска:

###  Интерактивный режим
Вы можете запускать скрипт с пошаговым вводом параметров вручную:

```bash
python PocketMaster.py
```

###  Режим с конфигурационным файлом
Для автоматического запуска с заранее заданными параметрами используйте YAML-файл конфигурации:

```bash
python PocketMaster.py --config config.yaml
```
📌 Пример конфигурационного файла см. в папке examples/config.yaml.

Дополнительно, при запуске в интерактивном режиме в конце работы создаётся файл run_config.txt, в который автоматически сохраняются все введенные
пользователем параметры. Это позволяет в дальнейшем выполнять запуск с теми же настройками в автоматическом режиме, без необходимости повторного ввода исходных данных.

  
## 📋 Как использовать

Запустите скрипт **PocketMaster.py** в интерактивном режиме командой, приведённой выше,  
чтобы перейти в интерактивный режим и задать параметры вручную.


🧭 Следуйте интерактивным подсказкам:

🔹 **Выбери режим работы:**  
  1 –  Использовать локальную папку с PDB файлами  
  2 –  На основе UniProt ID скачать все соответствующие PDB структуры  
  3 –  На основе PDB ID определить UniProt ID и скачать все соответствующие PDB структуры  

🔹 **Провести предварительную очистку структур?**    
  1 –  Да  
  2 –  Нет  

🔹 **Выберите подходящие опции предварительной обработки структур**    
  1 –  Удалить воду (solvent)  
  2 –  Удалить ионы ( Cu, CL,  и др.)  
  3 –  Удалить сульфаты и фосфаты (SO4, PO4, и др.)  
  4 –  Удалить буферные компоненты (TRS, MES, HEP, и др.)  
  5 –  Удалить криопротектанты (GOL, EDO, MPD, и др.)  
  6 –  Удалить восстановители (DTT, BME, TCEP)  
  7 –  Удалить всё воду, ионы, буферы, криопротектанты, фосфаты, восстановители   
  8 –  Удалить модифицированные аминокислотные остатки (CSO, MSE, SEP, TPO, PTR и др.)   
  9 –  Удалить всё, кроме белка (оставить только полимерную цепь   
  10 – Удалить альтернативные конформации (altloc)  
  11 – Удалить анизотропные параметры (ANISOU)  
  12 – Удалить атомы водорода (H)  
  13 – Сохранить обработанные структуры в определённой папке  
  14 – Не очищать / Завершить выбор  
  
🔹 **Выберите референсную структуру для выравнивания**  
  – введя её порядковый номер из представленного списка: 

🔹 **Выбор способа определения участка для выравнивания**  
  1 – На реф. структуре по заданному ID остатка и радиусу (Å), затем он ищется и выравнивается во всех структурах  
  2 – После предварительного выравнивания всех структур между собой, для каждой структуры вокруг выбранного реф. остатка  
  3 – Для каждой структуры вокруг её HET-групп в пределах заданного радиуса (Å)  
  4 – По введённому пользователем списку остатков, затем ищется во всех структурах  
  5 – На реф. структуре по указанному идентификатору цепи, затем ищется и выравнивается во всех структурах 

***
<p align="left">Поддерживаемые методы определения участков для выравнивания и их детальная интерпретация.</p>
<p align="left">
  <img src="images/graphic.png" height="400" />
</p>

***

🔹 **Выберите режим сравнения:**  
  1 – «все со всеми» (all vs all)  
  2 – «все с референсом» (all vs ref)  

🔹 **Укажите метод расчёта RMSD:**  
  1 – `align`строгая RMSD-оценка, используется для точного выравнивания, исключая атомы, не имеющие соответствий.  
  2 – `cealign` Структурное выравнивание на основе геометрии, эффективно даже при низком сходстве.  
  3 – `super` гибкий вариант, допускает несовпадения и автоматически подбирает соответствующие атомы.  
  4 – `rms` rms: точный RMSD, работает только при полном соответствии атомов (без несовпадений).  
  5 – `rms_cur` rms_cur: упрощённый и более быстрый расчёт RMSD, также требует полного соответствия атомов.  

🔹 **Выбери метод иерархической кластеризации:**  
  1 – ward (минимизация внутрикластерной дисперсии, требует евклидово расстояние)  
  2 – single (минимальное расстояние между кластерами)  
  3 – complete (максимальное расстояние между кластерами)  
  4 – average (среднее расстояние между кластерами (UPGMA))  
  5 – centroid (расстояние между центрами масс кластеров)  
  6 – median (медианное расстояние между кластерами)  
  7 – weighted (взвешенное среднее расстояние (WPGMA))  

<h3 align="left">📊 Сравнение методов кластеризации</h3>

<p align="left">

| Метод      | Компактные кластеры | Вытянутые кластеры | Чувствительность к выбросам |
|-------------|--------------------|--------------------|------------------------------|
| Ward        | ✅ Отлично          | ❌ Плохо           | ⚠️ Средняя                   |
| Single      | ❌ Плохо            | ✅ Отлично         | ⚠️ Высокая                   |
| Complete    | ✅ Хорошо           | ❌ Плохо           | ✅ Устойчив                  |
| Average     | ✅ Универсально     | ✅ Неплохо         | ⚠️ Средняя                   |
| Centroid    | ⚠️ Нестабильно      | ⚠️ Нестабильно     | ⚠️ Нестабильно               |
| Median      | ⚠️ Нестабильно      | ⚠️ Нестабильно     | ⚠️ Нестабильно               |
| Weighted    | ✅ Хорошо           | ✅ Хорошо          | ⚠️ Средняя                   |

</p>
<br>

🔹 **Параметры кластеризации: количество или порог**  
  1 – Задать количество кластеров (maxclust)  
  2 – Задать порог расстояния (distance)  
  3 – Автоматический: Порог на основе 70% процента максимального расстояния слиянияя  
  4 – Автоматический: Метод Изгиба или Локтя (Elbow Method)  

🔹 **Выберите тип RMSD:**  
  1 – по всем атомам  
  2 – только по Cα-атомам

## После выполнения

Скрипт создаёт папку `output` в указанной директории, содержащую:

1. `aligned_structures/`: Выровненные PDB файлы (например, `structure_aligned_to_ref.pdb`).
2. `rmsd_all_atoms.csv`: RMSD для всех атомов кармана.
3. `rmsd_calpha.csv`: RMSD только для Cα атомов.
4. `info.txt`: Информация о референсной структуре, методе задания кармана,  выравниваниях и сравнении остатков карманов.
5. `rmsd_all_atoms_hist.png`: Гистограмма RMSD для всех атомов.
6. `rmsd_calpha_hist.png`: Гистограмма RMSD для Cα атомов.
7. `Результаты кластеризации`: cluster_assignments.csv.
8. `Тепловая карта`: rmsd_heatmap.png.
9. `Дендрограмма`: RMSD_dendrogram.png.

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


