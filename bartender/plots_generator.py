import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def generate_demux_plots(report_tsv_path, output_dir):
    """
    Читает demux_report.tsv и строит графики распределения статусов и образцов.
    """
    # Настройка стиля для красивой публикации
    sns.set_theme(style="whitegrid")
    plt.rcParams.update({'font.size': 11, 'axes.labelsize': 12, 'axes.titlesize': 14})
    
    # 1. Загрузка данных отчета
    if not os.path.exists(report_tsv_path):
        print(f"Ошибка: Файл отчета {report_tsv_path} не найден.")
        return
        
    df = pd.read_csv(report_tsv_path, sep='\t', index_col=False)
    
    # --- ГРАФИК 1: Распределение статусов (Donut Chart) ---
    status_counts = df['status'].value_counts()
    
    fig, ax = plt.subplots(figsize=(7, 6))
    colors = sns.color_palette("pastel")[0:len(status_counts)]
    
    wedges, texts, autotexts = ax.pie(
        status_counts, 
        labels=status_counts.index, 
        autopct='%1.1f%%',
        startangle=90, 
        colors=colors, 
        textprops=dict(color="black"),
        pctdistance=0.75
    )
    
    # Превращаем круговую диаграмму в "пончик"
    centre_circle = plt.Circle((0,0), 0.55, fc='white')
    fig.gca().add_artist(centre_circle)
    
    # Делаем подписи процентов более читаемыми внутри секторов
    plt.setp(autotexts, size=10, weight="bold")
    ax.set_title("Распределение ридов по статусам фильтрации", pad=20)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "demux_status_distribution.png"), dpi=300)
    plt.close()

    # --- ГРАФИК 2: Распределение ридов по SampleID (Bar Chart) ---
    # Фильтруем только успешно определенные риды
    success_df = df[df['status'] == 'Fully_Identified']
    
    if success_df.empty:
        print("Предупреждение: Нет успешно определенных ридов для построения графика образцов.")
        return
        
    sample_counts = success_df['sample_id'].value_counts().reset_index()
    sample_counts.columns = ['SampleID', 'Read Count']
    
    # Ограничиваем топ-20 образцов, чтобы график не сливался, если образцов тысячи
    top_samples = sample_counts.head(20)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.barplot(
        data=top_samples, 
        x='Read Count', 
        y='SampleID', 
        hue='SampleID', 
        palette="viridis", 
        legend=False,
        ax=ax
    )
    
    ax.set_title("Топ-20 образцов по количеству ридов", pad=15)
    ax.set_xlabel("Количество прочтений (ридов)")
    ax.set_ylabel("Идентификатор образца (Sample ID)")
    
    # Добавляем текстовые значения количества ридов на концы столбцов
    for i, v in enumerate(top_samples['Read Count']):
        ax.text(v + (v * 0.01), i + .15, f" {v:,}", color='black', va='center', fontsize=9)
        
    plt.savefig(os.path.join(output_dir, "demux_samples_yield.png"), dpi=300)
    plt.close()
    
    print(f"-> Графики успешно сохранены в папку: '{output_dir}'")
