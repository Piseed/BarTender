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
        print(f"Error: Report file {report_tsv_path} not found.")
        return
        
    df = pd.read_csv(report_tsv_path, sep='\t', index_col=False)
    
    # --- ГРАФИК 1: Распределение статусов (Donut Chart) ---
    status_counts = df['status'].value_counts()
    
    fig, ax = plt.subplots(figsize=(7, 6))
    colors = sns.color_palette("cividis")[0:len(status_counts)]
    
    wedges, texts, autotexts = ax.pie(
        status_counts, 
        labels=status_counts.index, 
        autopct='%1.1f%%',
        rotatelabels=True,
        startangle=90, 
        colors=colors, 
        textprops={"color":"black", "fontsize":9},
        pctdistance=0.75
    )
    for text in texts:
      text.set_rotation(30)
    # Превращаем круговую диаграмму в "пончик"
    centre_circle = plt.Circle((0,0), 0.55, fc='white')
    fig.gca().add_artist(centre_circle)
    
    # Делаем подписи процентов более читаемыми внутри секторов
    plt.setp(autotexts, size=10, weight="bold")
    ax.set_title("Demux status distribution", pad=20)
    
    plt.savefig(os.path.join(output_dir, "demux_status_distribution.png"), dpi=300, bbox_inches='tight')
    plt.close()

    # --- ГРАФИК 2: Распределение ридов по SampleID (Bar Chart) ---
    # Фильтруем только успешно определенные риды
    success_df = df[df['status'] == 'Fully_Identified']
    
    if success_df.empty:
        print("No fully identified reads. The bar plot will not be drawn.")
        return
        
    sample_counts = success_df['sample_id'].value_counts().reset_index()
    sample_counts.columns = ['SampleID', 'Read Count']
    sample_counts['SampleID'] = sample_counts['SampleID'].astype(int)
    # Ограничиваем топ-20 образцов, чтобы график не сливался, если образцов тысячи
    top_samples = sample_counts.head(20)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.barplot(
        data=top_samples, 
        x='SampleID', 
        y='Read Count', 
        hue='SampleID', 
        palette="cividis", 
        legend=False,
        ax=ax
    )
    
    ax.set_title("Тоp-20 samples by read count", pad=15)
    ax.set_xlabel("Read count")
    ax.set_ylabel("Sample ID")
    
    # Добавляем текстовые значения количества ридов на концы столбцов
        
    plt.savefig(os.path.join(output_dir, "demux_samples_yield.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"-> Plots saved to directory: '{output_dir}'")
