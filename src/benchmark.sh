#!/bin/bash
# Benchmark script for k-mer cache optimization

echo "=== K-mer Cache Optimization Benchmark ==="
echo ""

# Run multiple times and average
RUNS=10

echo "Running $RUNS iterations..."
echo ""

total_time=0
for i in $(seq 1 $RUNS); do
    start=$(python3 -c 'import time; print(time.time())')
    ./hmm --sequence seq \
        --don_emission ../models/don.pwm \
        --acc_emission ../models/acc.pwm \
        --exon_emission ../models/exon.mm \
        --intron_emission ../models/intron.mm \
        --ped_exon ../models/exon.len \
        --ped_intron ../models/intron.len \
        --n_isoforms 1000 \
        --mtry 0.5 \
        --stovit > /dev/null 2>&1
    end=$(python3 -c 'import time; print(time.time())')
    elapsed=$(python3 -c "print(f'{($end - $start):.4f}')")
    echo "Run $i: ${elapsed}s"
    total_time=$(python3 -c "print($total_time + $elapsed)")
done

avg_time=$(python3 -c "print(f'{$total_time / $RUNS:.4f}')")
echo ""
echo "Average time: ${avg_time}s"
