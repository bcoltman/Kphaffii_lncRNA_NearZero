set -ue

# Get the bam file from the command line
DATA=$1
NAME=$2
OUTPUT_DIR=$3

mkdir -p $OUTPUT_DIR 


TMPDIR=$(mktemp -d)

# Forward strand.
#
# 1. alignments of the second in pair if they map to the forward strand
# 2. alignments of the first in pair if they map to the reverse  strand
#
samtools view -b -f 128 -F 16 --threads 6 $DATA > $TMPDIR/${NAME}_fwd1.bam
samtools index -@ 6 $TMPDIR/${NAME}_fwd1.bam

samtools view -b -f 80 --threads 6 $DATA > $TMPDIR/${NAME}_fwd2.bam
samtools index -@ 6 $TMPDIR/${NAME}_fwd2.bam

#
# Combine alignments that originate on the forward strand.
#
samtools merge --threads 6 -f $OUTPUT_DIR/${NAME}_fwd.bam $TMPDIR/${NAME}_fwd1.bam $TMPDIR/${NAME}_fwd2.bam
samtools index -@ 6 $OUTPUT_DIR/${NAME}_fwd.bam

# Reverse strand
#
# 1. alignments of the second in pair if they map to the reverse strand
# 2. alignments of the first in pair if they map to the forward strand
#
samtools view -b -f 144 --threads 6 $DATA > $TMPDIR/${NAME}_rev1.bam
samtools index -@ 6 $TMPDIR/${NAME}_rev1.bam

samtools view -b -f 64 -F 16 --threads 6 $DATA > $TMPDIR/${NAME}_rev2.bam
samtools index -@ 6 $TMPDIR/${NAME}_rev2.bam

#
# Combine alignments that originate on the reverse strand.
#
samtools merge --threads 6 -f $OUTPUT_DIR/${NAME}_rev.bam $TMPDIR/${NAME}_rev1.bam $TMPDIR/${NAME}_rev2.bam
samtools index -@ 6 $OUTPUT_DIR/${NAME}_rev.bam

rm -r $TMPDIR
