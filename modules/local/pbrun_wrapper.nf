def pbrunFunction() {
    '''
pbrun() {
    local rc=0 exe
    exe=$(type -P pbrun)
    "$exe" "$@" 2>&1 | tee -a pbrun.log || rc=$?
    if [ "$rc" -ne 0 ] && grep -q 'SIGKILL' pbrun.log; then
        echo "pbrun was killed for host memory (SIGKILL); exiting 137 so the task is retried with more memory" >&2
        return 137
    fi
    return "$rc"
}
'''
}

def dedupeHdFunction() {
    '''
dedupe_hd() {
    local f=$1 n
    n=$(samtools view -H "$f" | grep -c '^@HD' || true)
    [ "$n" -gt 1 ] || return 0
    if ! samtools view -H "$f" | grep -q '^@HD.*SO:coordinate'; then
        echo "$f has $n @HD lines and none is SO:coordinate; not rewriting the header" >&2
        return 1
    fi
    echo "$f has $n @HD lines; keeping only the SO:coordinate line" >&2
    local prog='/^@HD/ { if (done || $0 !~ /SO:coordinate/) next; done=1 } { print }'
    samtools reheader --no-PG -c "awk '$prog'" "$f" > "$f.reheadered"
    mv "$f.reheadered" "$f"
    rm -f "$f.bai" "$f.crai" "${f%.*}.bai"
}
'''
}
