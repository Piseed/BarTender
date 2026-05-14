INPUT_FASTA="${1:-barcodes.fasta}"
LOG_FILE="pipeline_manifest.log"

echo "=== [$(date)] Запуск пайплайна ===" > "$LOG_FILE"
echo "Исходный файл: $INPUT_FASTA" >> "$LOG_FILE"

touch tmp_F.fasta tmp_R.fasta tmp_discard.fasta
# --- Шаг 1: Разделение на группы F и R ---
awk '
    /^>/ { f = (substr($1,2,1) == "F") ? "tmp_F.fasta" : ((substr($1,2,1) == "R") ? "tmp_R.fasta" : "tmp_discard.fasta") }
    { print $0 >> f }
' "$INPUT_FASTA"

# --- Функция поиска адаптера с 3-конца ---
find_3_adapter() {
    if [ ! -s "$1" ] || [ $(awk 'NR%2==0' "$1" | wc -l) -eq 0 ]; then
        echo ""
        return
    fi
    awk 'NR%2==0' "$1" | awk '
        NR==1 { res = $0; next }
        {
            while (substr($0, length($0)-length(res)+1) != res) {
                res = substr(res, 2)
            }
            if (res == "") exit
        }
        END { print res }
    '
}

ADAPTER_F=$(find_3_adapter "tmp_F.fasta")
ADAPTER_R=$(find_3_adapter "tmp_R.fasta")

echo "Найден адаптер для группы F: $ADAPTER_F" >> "$LOG_FILE"
echo "Найден адаптер для группы R: $ADAPTER_R" >> "$LOG_FILE"

# --- Функция обработки группы (Только RAW и REVCOMP) ---
process_group() {
    local group_name="$1"    # "F" или "R"
    local tmp_file="$2"      # "tmp_F.fasta"
    local adapter="$3"       # Последовательность адаптера
    local ad_len=${#adapter}

    echo "Обработка группы $group_name..."
    if [ ! -s "$tmp_file" ]; then
        return
    fi
    # 1. Извлекаем сырую коллекцию (удаляем адаптер с 3-конца)
    awk -v len="$ad_len" '
        NR%2!=0 { print $0 }
        NR%2==0 { print substr($0, 1, length($0) - len) }
    ' "$tmp_file" > "collection_${group_name}_raw.fasta"

    # Генерация обратного комплемента для адаптера
    local ad_revcomp=""
    if [ -n "$adapter" ]; then
        ad_revcomp=$(echo "$adapter" | tr ATCGatcg TAGCtagc | rev)
    fi

    # Запись адаптеров в манифест
    echo "  - Adapter ${group_name}_raw: $adapter" >> "$LOG_FILE"
    echo "  - Adapter ${group_name}_revcomp: $ad_revcomp" >> "$LOG_FILE"

    # 2. Генерируем коллекцию обратно-комплементарных баркодов
    awk '
        NR%2!=0 { print $0 "_revcomp" }
        NR%2==0 {
            seq = $0
            cmd = "echo " seq " | tr ATCGatcg TAGCtagc"
            cmd | getline seq
            close(cmd)
            
            split(seq, chars, ""); res = "";
            for (i=length(seq); i>0; i--) res = res chars[i];
            print res
        }
    ' "collection_${group_name}_raw.fasta" > "collection_${group_name}_revcomp.fasta"
}

# Запуск генерации для обеих групп
process_group "F" "tmp_F.fasta" "$ADAPTER_F"
process_group "R" "tmp_R.fasta" "$ADAPTER_R"

# --- Автоматический экспорт конфигурации адаптеров в JSON ---
echo "{\"adapterf\": \"$ADAPTER_F\", \"adapterr\": \"$ADAPTER_R\"}" > info_adapters.json

# Логирование созданных файлов
echo -e "\nСозданные файлы коллекций баркодов:" >> "$LOG_FILE"
ls -1 collection_*.fasta >> "$LOG_FILE"

# Очистка
rm -f tmp_F.fasta tmp_R.fasta tmp_discard.fasta

echo "=== Пайплайн успешно завершен ==="
echo "Все действия зафиксированы в: $LOG_FILE"