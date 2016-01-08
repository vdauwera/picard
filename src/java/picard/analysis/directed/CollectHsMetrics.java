/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package picard.analysis.directed;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.StringUtil;
import picard.analysis.MetricAccumulationLevel;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Metrics;

import java.io.File;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * Collects a set of HS metrics from a sam or bam file.  See HsMetricsCollector and CollectTargetedMetrics for more details.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        usage = CollectHsMetrics.USAGE_SUMMARY + CollectHsMetrics.USAGE_DETAILS,
        usageShort = CollectHsMetrics.USAGE_SUMMARY,
        programGroup = Metrics.class

)
public class CollectHsMetrics extends CollectTargetedMetrics<HsMetrics, HsMetricCollector> {

    static final String USAGE_SUMMARY = "Collects hybrid-selection (HS) specific metrics for a SAM or BAM file.  ";
    static final String USAGE_DETAILS = "These metrics enable users to determine the efficacy of HS experiments.  " +
            "In addition, if a reference sequence is provided, this program will calculate AT/GC dropout metrics and enable users to invoke the" +
            " PER_TARGET_COVERAGE option to output both GC and mean coverage information for every target. <br /><br />" +
            "" +
            "Both CollectTargetedPCRMetrics and CalculateHybridSelection metrics share virtually identical program structures except " +
            "for the name of their targeting mechanisms (e.g. bait set or amplicon set).  Both depend on the TargetMetricsCollector " +
            "and collect metrics on all reads in the INPUT SAM/BAM file. <br /><br /> " +
            "" +
            "Hybrid selection (HS) enables targeted sequencing analysis via the capture of specified genomic DNA " +
            "sequences <a href=\"http://www.nature.com/nbt/journal/v27/n2/abs/nbt.1523.html\"><strong>(doi:10.1038/nbt.1523)</strong></a>.  " +
            "It is commonly used to characterize exon sequences from genomic DNA " +
            "or filter out bacterial DNA sequences from clinical samples.  <br /><br /> " +
            "" +
            "The technique involves the capture of unique regions of genomic DNA (targets) using synthetic RNA baits." +
            "   The baits are synthesized with biotinylated nucleotides to facilitate capture of bait:target hybrids on" +
            " streptavidin beads.  The captured target sequences are amplified, sequenced, and processed for variant calling. <br /><br />" +
            "" +
            "This tool requires an aligned SAM or BAM file as well as bait and target interval_list files.  " +

            "For information on interval lists, please see the documentation for IntervalListTools" +
            "<a href=\"http://broadinstitute.github.io/picard/command-line-overview.html#IntervalListTools\"><strong> here</strong>.  </a>" +
            "<p>CollectHsMetrics provides multiple outputs that are described in detail at the following link: " +
            "<a href=\"http://broadinstitute.github.io/picard/picard-metric-definitions.html\"><strong>metrics definitions</strong></a>.</p>  " +
            ""+
            "<h4>Usage Example:</h4>"+
            "<pre>" +
            "java -jar picard.jar CollectHsMetrics \\<br />" +
            "      I=input.bam \\<br />" +
            "      O=hs_metrics.txt \\<br />" +
            "      R=reference_sequence.fasta \\<br />" +
            "      BAIT_INTERVALS=bait.interval_list \\<br />" +
            "      TARGET_INTERVALS=target.interval_list" +
            "</pre> "+
            "<hr />";


    @Option(shortName = "BI", doc = "An interval list file that contains the locations of the baits used.", minElements=1)
    public List<File> BAIT_INTERVALS;

    @Option(shortName = "N", doc = "Bait set name. If not provided it is inferred from the filename of the bait intervals.", optional = true)
    public String BAIT_SET_NAME;

    @Option(shortName = "MQ", doc = "Minimum mapping quality for a read to contribute coverage.", overridable = true)
    public int MINIMUM_MAPPING_QUALITY = 20;

    @Option(shortName = "Q", doc = "Minimum base quality for a base to contribute coverage.", overridable = true)
    public int MINIMUM_BASE_QUALITY = 20;

    @Option(doc = "True if we are to clip overlapping reads, false otherwise.", optional=true, overridable = true)
    public boolean CLIP_OVERLAPPING_READS = true;

    @Override
    protected IntervalList getProbeIntervals() {
        for (final File file : BAIT_INTERVALS) IOUtil.assertFileIsReadable(file);
        return IntervalList.fromFiles(BAIT_INTERVALS);
    }

    @Override
    protected String getProbeSetName() {
        if (BAIT_SET_NAME != null) {
            return BAIT_SET_NAME;
        } else {
            final SortedSet<String> baitSetNames = new TreeSet<String>();
            for (final File file : BAIT_INTERVALS) {
                baitSetNames.add(CollectTargetedMetrics.renderProbeNameFromFile(file));
            }
            return StringUtil.join(".", baitSetNames);
        }
    }

    /** Stock main method. */
    public static void main(final String[] argv) {
        System.exit(new CalculateHsMetrics().instanceMain(argv));
    }

    @Override
    protected HsMetricCollector makeCollector(final Set<MetricAccumulationLevel> accumulationLevels,
                                              final List<SAMReadGroupRecord> samRgRecords,
                                              final ReferenceSequenceFile refFile,
                                              final File perTargetCoverage,
                                              final File perBaseCoverage,
                                              final IntervalList targetIntervals,
                                              final IntervalList probeIntervals,
                                              final String probeSetName,
                                              final int nearProbeDistance) {
        return new HsMetricCollector(accumulationLevels, samRgRecords, refFile, perTargetCoverage, perBaseCoverage, targetIntervals, probeIntervals, probeSetName, nearProbeDistance,
                MINIMUM_MAPPING_QUALITY, MINIMUM_BASE_QUALITY, CLIP_OVERLAPPING_READS, true, COVERAGE_CAP, SAMPLE_SIZE);
    }
}