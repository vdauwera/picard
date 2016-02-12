/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
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
package picard.analysis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.OverlapDetector;
import picard.PicardException;
import picard.analysis.directed.RnaSeqMetricsCollector;
import picard.annotation.Gene;
import picard.annotation.GeneAnnotationReader;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Metrics;
import picard.util.RExecutor;

import java.io.File;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

@CommandLineProgramProperties(
        usage = CollectRnaSeqMetrics.USAGE_SUMMARY + CollectRnaSeqMetrics.USAGE_DETAILS,
        usageShort = CollectRnaSeqMetrics.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class CollectRnaSeqMetrics extends SinglePassSamProgram {

    static final String USAGE_SUMMARY = "Produces RNA alignment metrics for a SAM or BAM file.  ";
    static final String USAGE_DETAILS = "<p>Program takes aligned RNA transcript sequences and determines the functional distribution of bases " +
            "within the transcripts of a BAM file including the numbers and percentages of bases that are located within " +
            "the untranslated region (UTR), introns, intergenic (sequences between " +
            "discrete genes), as well peptide-coding sequences.  This tool also calculates quality metrics on the reads and the bases within" +
            " the reads including the numbers of bases that pass the sequencing vendor's quality filters (PF_BASES).  " +
            "For additional information on the PF metric, please see the GATK dictionary entry on" +
            "<a href=\"http://gatkforums.broadinstitute.org/gatk/discussion/6329/pf-reads-illumina-chastity-filter\"><strong> (PF Reads/Illumina Chastity Filters)</strong>.  </a></p>" +
            "" +
            "Other metrics include the median coverage (depth), 5&apos;/3&apos; biases, and the numbers of reads with the correct/incorrect strand " +
            "designation.  The 5&apos;/3&apos; biases represent the bias introduced by having a" +
            " greater representation at one end of a transcript over another.  This bias is most often introduced during library construction resulting from " +
            "e.g. incomplete activity of a reverse transcriptase.  For additional information, please see the following" +
            "" +
            "<a href=\"http://www.biostars.org/p/102812/#102815\"><strong> link</strong>.  </a>" +
            "" +
            "<p>Users must have a functional BAM file containing RNA transcripts generated from an aligner such as <a href=\"http://github.com/alexdobin/STAR\"><strong> STAR</strong></a>" +
            " or <a href=\"http://github.com/infphilo/tophat\"><strong> TopHat</strong></a>.  Prior to input into CollectRnaSeqMetrics, BAM files should be validated using Picard's" +
            "" +
            "<a href=\"http://broadinstitute.github.io/picard/command-line-overview.html#ValidateSamFile\"><strong> ValidateSamFile</strong></a> tool."+
            "" +
            "" +
            "" +
            "<p>CollectRnaSeqMetrics also requires REF_FLAT file," +
            " a tab-delimited file containing information about the location of " +
            "RNA transcripts, exon start and stop sites, etc.  Please visit the following link for documentation on " +
            " <a href=\"http://genome.ucsc.edu/FAQ/FAQformat.html\"> <strong>REF_FLAT</strong> </a> files and search the page for the phrase \"refFlat\".  " +
            "" +
            " To obtain genome-specific downloadable REF_FLAT files, please visit the link " +
            "<a href=\"http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/\"> <strong>here</strong></a>.</p>" +
            "" +
            "<p>If provided with an interval file identifying ribosomal sequences, it calculates numbers/percentages of ribosomal bases.  " +
            "The following site contains a ribosomal sequence interval list for the " +
            "<a href=\"http://gist.github.com/slowkow/b11c28796508f03cdf4b\"><strong>Hg19</strong></a> genomic build.  " +
            "Note, you can copy and paste the contents from this file into a text editor to create your own ribosomal sequence interval file.  In addition, if you are not " +
            "using Hg19, you can use Picard's " +
            "" +
            " <a href=\"http://broadinstitute.github.io/picard/command-line-overview.html#LiftOverIntervalList\"><strong> LiftOverIntervalList</strong></a> tool to lift the file over to another genomic build. </p>" +
            "" +
            "<pre>" +
            "<h4>Usage example:</h4>" +
            "java -jar picard.jar CollectRnaSeqMetrics \\<br />" +
            "      I=input.bam \\<br />" +
            "      O=output.RNA_Metrics \\<br />" +
            "      REF_FLAT=ref_flat.txt \\<br />" +
            "      STRAND=SECOND_READ_TRANSCRIPTION_STRAND \\<br />" +
            "      RIBOSOMAL_INTERVALS=ribosomal.interval_list <br />" +
            "</pre>" +


            "<p>Users can invoke the &quotIGNORE_SEQUENCE&quot option to specify reads not be included in the analysis.</p>" +
            "" +
            "<p>For a complete description of the output metrics, please see " +
            " <a href=\"http://broadinstitute.github.io/picard/picard-metric-definitions.html\"><strong>metrics definitions</strong></a>.</p>  " +
            ""+
    "<hr />";

    private static final Log LOG = Log.getInstance(CollectRnaSeqMetrics.class);

    @Option(doc="Gene annotations in refFlat form.  Format described here: http://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#RefFlat")
    public File REF_FLAT;

    @Option(doc="Location of rRNA sequences in genome, in interval_list format.  " +
            "If not specified no bases will be identified as being ribosomal.  " +
            "Format described here: http://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/util/IntervalList.html", optional = true)
    public File RIBOSOMAL_INTERVALS;

    @Option(shortName = "STRAND", doc="For strand-specific library prep. " +
            "For unpaired reads, use FIRST_READ_TRANSCRIPTION_STRAND if the reads are expected to be on the transcription strand.")
    public RnaSeqMetricsCollector.StrandSpecificity STRAND_SPECIFICITY;

    @Option(doc="When calculating coverage based values (e.g. CV of coverage) only use transcripts of this length or greater.")
    public int MINIMUM_LENGTH = 500;

    @Option(doc="The PDF file to write out a plot of normalized position vs. coverage.", shortName="CHART", optional = true)
    public File CHART_OUTPUT;

    @Option(doc="If a read maps to a sequence specified with this option, all the bases in the read are counted as ignored bases.  " +
    "These reads are not counted as ")
    public Set<String> IGNORE_SEQUENCE = new HashSet<String>();

    @Option(doc="This percentage of the length of a fragment must overlap one of the ribosomal intervals for a read or read pair by this must in order to be considered rRNA.")
    public double RRNA_FRAGMENT_PERCENTAGE = 0.8;

    @Option(shortName="LEVEL", doc="The level(s) at which to accumulate metrics.  ")
    public Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = CollectionUtil.makeSet(MetricAccumulationLevel.ALL_READS);

    private RnaSeqMetricsCollector collector;

    /**
     * A subtitle for the plot, usually corresponding to a library.
     */
    private String plotSubtitle = "";

    /** Required main method implementation. */
    public static void main(final String[] argv) {
        new CollectRnaSeqMetrics().instanceMainWithExit(argv);
    }

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {

        if (CHART_OUTPUT != null) IOUtil.assertFileIsWritable(CHART_OUTPUT);

        final OverlapDetector<Gene> geneOverlapDetector = GeneAnnotationReader.loadRefFlat(REF_FLAT, header.getSequenceDictionary());
        LOG.info("Loaded " + geneOverlapDetector.getAll().size() + " genes.");

        final Long ribosomalBasesInitialValue = RIBOSOMAL_INTERVALS != null ? 0L : null;
        final OverlapDetector<Interval> ribosomalSequenceOverlapDetector = RnaSeqMetricsCollector.makeOverlapDetector(samFile, header, RIBOSOMAL_INTERVALS);

        final HashSet<Integer> ignoredSequenceIndices = RnaSeqMetricsCollector.makeIgnoredSequenceIndicesSet(header, IGNORE_SEQUENCE);

        collector = new RnaSeqMetricsCollector(METRIC_ACCUMULATION_LEVEL, header.getReadGroups(), ribosomalBasesInitialValue,
                geneOverlapDetector, ribosomalSequenceOverlapDetector, ignoredSequenceIndices, MINIMUM_LENGTH, STRAND_SPECIFICITY, RRNA_FRAGMENT_PERCENTAGE,
                true);

        // If we're working with a single library, assign that library's name as a suffix to the plot title
        final List<SAMReadGroupRecord> readGroups = header.getReadGroups();
        if (readGroups.size() == 1) {
            this.plotSubtitle = readGroups.get(0).getLibrary();
            if (null == this.plotSubtitle) this.plotSubtitle = "";
        }
    }

    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence refSeq) {
        collector.acceptRecord(rec, refSeq);
    }

    @Override
    protected void finish() {
        collector.finish();

        final MetricsFile<RnaSeqMetrics, Integer> file = getMetricsFile();
        collector.addAllLevelsToFile(file);
        file.write(OUTPUT);

        boolean atLeastOneHistogram = false;
        for (final Histogram<Integer> histo : file.getAllHistograms()) {
            atLeastOneHistogram = atLeastOneHistogram || !histo.isEmpty();
        }
        // Generate the coverage by position plot
        if (CHART_OUTPUT != null && atLeastOneHistogram) {
            final int rResult = RExecutor.executeFromClasspath("picard/analysis/rnaSeqCoverage.R",
                                                               OUTPUT.getAbsolutePath(),
                                                               CHART_OUTPUT.getAbsolutePath(),
                                                               INPUT.getName(),
                                                               this.plotSubtitle);

            if (rResult != 0) {
                throw new PicardException("Problem invoking R to generate plot.");
            }
        }
    }

}
