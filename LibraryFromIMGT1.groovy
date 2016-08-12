@Grab(group = 'com.milaboratory', module = 'repseqio', version = '1.0.1-SNAPSHOT')
import cc.redberry.pipe.CUtils
import com.milaboratory.core.Range
import com.milaboratory.core.io.sequence.fasta.FastaReader
import com.milaboratory.core.io.sequence.fasta.FastaRecord
import com.milaboratory.core.sequence.NucleotideSequence
import com.milaboratory.util.GlobalObjectMappers
import io.repseq.core.BaseSequence
import io.repseq.core.Chains
import io.repseq.core.GeneFeature
import io.repseq.core.ReferencePoint
import io.repseq.dto.VDJCDataUtils
import io.repseq.dto.VDJCGeneData
import io.repseq.dto.VDJCLibraryData

if (args.size() < 3) {
    println "Usage: script.groovy input.fasta output.json taxonId common_species_name..."
    return
}

Map<String, VDJCGeneData> genes = new HashMap<>();

List<String> chains = ["TRA", "TRB", "TRG", "TRD", "IGH", "IGL", "IGK"];

def reader = new FastaReader<>(args[0], NucleotideSequence.ALPHABET)
for (FastaRecord<NucleotideSequence> record in CUtils.it(reader)) {
    def split = record.description.split("\\|")
    def rangeSplit = split[2].split("\\.\\.")
    int from = Integer.parseInt(rangeSplit[0]);
    int to = Integer.parseInt(rangeSplit[1]);
    def range = from < to ? new Range(from - 1, to) : new Range(from, to - 1)
    def geneName = split[4] + "*00"
    def accession = split[1]
    def featureS = split[5]
    def feature = null
    switch (featureS) {
        case "L-PART1":
            feature = GeneFeature.L1
            break
        case "V-INTRON":
            feature = GeneFeature.VIntron
            break
        case "V-EXON":
            feature = GeneFeature.VExon2
            break
        case "J-REGION":
            feature = GeneFeature.JRegion
            break
    }

    if (feature == null)
        continue

    // Assert simple (not complex) feature
    assert feature.size() == 1

    def gene = genes.get(geneName)
    if (gene == null) {
        def geneType = feature.geneType
        genes.put(geneName, gene = new VDJCGeneData(new BaseSequence("nuccore://" + accession), geneName, geneType, true,
                new Chains(chains.findAll { c -> geneName.indexOf(c[-1] + geneType.letter) >= 0 }.toSet()),
                new TreeMap<ReferencePoint, Long>()))
    }

    def anchorPoints = gene.anchorPoints

    if (anchorPoints.containsKey(feature.getReferenceRange(0).begin)) {
        assert anchorPoints.get(feature.getReferenceRange(0).begin).intValue() == range.from
//        println "OK"
    } else
        anchorPoints.put(feature.getReferenceRange(0).begin, new Long(range.from))

    if (anchorPoints.containsKey(feature.getReferenceRange(0).end)) {
        assert anchorPoints.get(feature.getReferenceRange(0).end).intValue() == range.to
//        println "OK"
    } else
        anchorPoints.put(feature.getReferenceRange(0).end, new Long(range.to))
}

VDJCLibraryData library = new VDJCLibraryData(Long.decode(args[2]), new ArrayList<String>(args[3..-1]),
        new ArrayList<VDJCGeneData>(genes.values()), Collections.EMPTY_LIST);

VDJCDataUtils.sort(library)

GlobalObjectMappers.PRETTY.writeValue(new File(args[1]), [library]);

//def f = VDJCLibraryRegistry.getDefault().getLibrary("default", "hs").get("TRBV12-2*00").getFeature(GeneFeature.VRegion)
//println f