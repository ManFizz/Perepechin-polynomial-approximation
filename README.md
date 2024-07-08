
## Исследование оперативности и адекватности полиномиального приближения функционально заданных данных с помощью SPMD и SIMD параллелизма
### by Перепечин Владимир Б9120-09.03.04
## Описание
Данный проект позволяет собирать результаты с использованием C++ и построение графиков с использованием Python и Jupyter Notebooks.

## Требования
Для успешного запуска проекта, необходимо установить следующие зависимости:

### Для запуска сбора результатов:
- C++20
- CMake версии 3.10 или выше
- MinGW64
- Библиотеку AVX-512 by Нельбасов Денис

### Для запуска построения графиков:
- Python
- Jupyter Notebooks

## Установка и настройка

### 1. Установка CMake и MinGW64
#### Windows:
1. **Установка CMake**:
    - Скачайте CMake с [официального сайта](https://cmake.org/download/).
    - Установите CMake, следуя инструкциям установщика.

2. **Установка MinGW64**:
    - Скачайте MinGW64 с [SourceForge](https://sourceforge.net/projects/mingw-w64/).
    - Установите MinGW64, следуя инструкциям установщика.
    - Добавьте путь к MinGW64 в переменную окружения PATH. Обычно это `C:\Program Files\mingw-w64\x86_64-<версия>\bin`.

3. **Установка AVX-512 библиотеки**:
   - Скачайте библиотеку с [GitHub](https://github.com/NellBaZZ/Diplom).
   - Создайте папку lib в директории этого проекта.
   - Разархивируйте библиотеку в созданную папку lib.

### 2. Установка Python и Jupyter Notebooks
#### Windows:
1. **Установка Python**:
    - Скачайте Python с [официального сайта](https://www.python.org/downloads/).
    - Установите Python, выбрав опцию "Add Python to PATH" во время установки.

2. **Установка Jupyter Notebooks**:
    - Откройте командную строку (cmd).
    - Установите Jupyter с помощью pip:
      ```sh
      pip install notebook
      ```

### 3. Сборка и запуск проекта

1. **Сборка C++ проекта**:
    - Откройте командную строку.
    - Перейдите в каталог проекта:
      ```sh
      cd /pathToProject/
      ```
    - Создайте директорию сборки и перейдите в нее:
      ```sh
      mkdir build
      cd build
      ```
    - Запустите CMake и сборку проекта:
      ```sh
      cmake ..
      cmake --build .
      ```

2. **Запуск Jupyter Notebooks**:
    - Откройте командную строку.
    - Перейдите в каталог с Jupyter Notebooks:
      ```sh
      cd path/to/your/notebooks
      ```
    - Запустите Jupyter Notebook:
      ```sh
      jupyter notebook
      ```
    - Откроется браузер с интерфейсом Jupyter, где вы сможете построить графики.

## Пример использования
### Сбор результатов
Запустите собранный C++ проект из директории сборки:
```sh
./cmake-build-*/VKR.exe
```

### Построение графиков
Запустите Jupyter Notebook и выполните все ячейки для построения графиков на основе собранных данных.