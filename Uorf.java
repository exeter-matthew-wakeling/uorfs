import htsjdk.samtools.reference.*;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Calculate the change in uOrfs in a 5-prime UTR of a gene.
 *
 * @author Matthew Wakeling
 */
public class Uorf implements Comparable<Uorf>
{
	/**
	 * Calculate the uorf consequence for a variant, according to the command-line arguments.
	 * <ul><li>A fasta file containing the reference genome.</li>
	 *     <li>The chromosome the variant is in.</li>
	 *     <li>The position of the variant.</li>
	 *     <li>The reference allele.</li>
	 *     <li>The alternate allele.</li>
	 *     <li>The strand direction for the gene transcript: 1 for forward, and -1 for reverse.</li>
	 * </ul>
	 * This is then followed by a list of 5-prime UTR exons. Each one has two arguments, which are the start and end positions (inclusive).
	 * <p>
	 * For instance:<br>
	 * java Uorf human_g1k_v37.fasta 19 633529 G GGCGCCGCCGCCGCCGCCGCC -1 633513 633568
	 *
	 * @param args the command-line arguments
	 */
	public static void main(String[] args) throws Exception {
		IndexedFastaSequenceFile reference = new IndexedFastaSequenceFile(new File(args[0]));
		String chr = args[1];
		int pos = Integer.parseInt(args[2]);
		String ref = args[3];
		String alt = args[4];
		boolean forward = Integer.parseInt(args[5]) > 0;
		List<FivePrimeUtrExon> exons = new ArrayList<FivePrimeUtrExon>();
		for (int i = 6; i < args.length; i += 2) {
			exons.add(new FivePrimeUtrExon(chr, Integer.parseInt(args[i]), Integer.parseInt(args[i + 1])));
		}
		FivePrimeUtr fivep = new FivePrimeUtr(forward, exons);
		UorfResult result = calculateUorfEffect(reference, fivep, chr, pos, ref, alt);
		System.out.println("Effect: " + result.getEffect());
		System.out.println("Start codon strength: " + result.getUorf().getStrengthString());
		System.out.println("Start codon distance: " + result.getUorf().getDistance());
		System.out.println("ORF finish distance: " + result.getUorf().getStopDistance());
		System.out.println("uORFs in reference: " + result.getRefUorfs());
		System.out.println("uORFs in alternate: " + result.getAltUorfs());
		System.out.println("Reference UTR bases, offset 0: " + result.getVisualisations()[0]);
		System.out.println("Alternate UTR bases, offset 0: " + result.getVisualisations()[1]);
		System.out.println("Reference UTR bases, offset 1: " + result.getVisualisations()[2]);
		System.out.println("Alternate UTR bases, offset 1: " + result.getVisualisations()[3]);
		System.out.println("Reference UTR bases, offset 2: " + result.getVisualisations()[4]);
		System.out.println("Alternate UTR bases, offset 2: " + result.getVisualisations()[5]);
	}

	/**
	 * Object representing the location of the 5-prime UTR, which can be spread over multiple exons.
	 */
	public static class FivePrimeUtr
	{
		private boolean forwardStrand;
		private List<FivePrimeUtrExon> exons;

		/**
		 * Creates a new 5-prime UTR.
		 *
		 * @param forwardStrand true if the gene is on the forward strand of the chromosome, false for the reverse strand
		 * @param exons a List of FivePrimeUtrExon objects describing (in exon number order) the regions that are part of the 5-prime UTR
		 */
		public FivePrimeUtr(boolean forwardStrand, List<FivePrimeUtrExon> exons) {
			this.forwardStrand = forwardStrand;
			this.exons = exons;
		}

		public boolean getForwardStrand() {
			return forwardStrand;
		}

		public List<FivePrimeUtrExon> getExons() {
			return exons;
		}
	}

	/**
	 * Object representing a section of a 5-prime UTR - that is, the UTR section of an exon.
	 */
	public static class FivePrimeUtrExon
	{
		private String chr;
		private int start, end;

		/**
		 * Creates an exon part of a 5-prime UTR.
		 *
		 * @param chr the chromosome the exon is in
		 * @param start the start position of the exon (inclusive)
		 * @param end the end position of the exon (inclusive), or the last position of the UTR before the start of the coding part of the gene
		 */
		public FivePrimeUtrExon(String chr, int start, int end) {
			this.chr = chr;
			this.start = start;
			this.end = end;
		}

		public String getChr() {
			return chr;
		}

		public int getStart() {
			return start;
		}

		public int getEnd() {
			return end;
		}
	}

	/**
	 * Calculate the uORFS in the 5-prime UTR in both reference and variant cases, then calculate any differences between the two.
	 *
	 * @param reference a htsjdk.samtools.reference.IndexedFastaSequenceFile object to allow the reference genome to be read
	 * @param fivePrimeUtr a FivePrimeUtr object describing where the UTR is
	 * @param chr the chromosome of the variant
	 * @param pos the position of the variant
	 * @param ref the reference allele in the area where the variant is
	 * @param alt the alternate allele in the area where the variant is
	 *
	 * @return a UorfResult object containing the results
	 */
	public static UorfResult calculateUorfEffect(IndexedFastaSequenceFile reference, FivePrimeUtr fivePrimeUtr, String chr, int pos, String ref, String alt) {
		if (fivePrimeUtr == null) {
			// Cannot create a result without a FivePrimeUtr
			return new UorfResult("", false, null, null, null, null);
		}
		FivePrimeUtrExon overlaps = null;
		// Find the exon that contains the variant
		for (FivePrimeUtrExon exon : fivePrimeUtr.getExons()) {
			if (exon.getChr().equals(chr) && (exon.getStart() <= pos) && (exon.getEnd() >= pos + ref.length() - 1)) {
				overlaps = exon;
			}
		}
		if (overlaps == null) {
			// Variant is not in the 5-prime UTR
			return new UorfResult("", false, null, null, null, null);
		}
		String refBases = "";
		String altBases = "";
		// Build a base string for the 5-prime UTR before and after the variant, after splicing
		if (fivePrimeUtr.getForwardStrand()) {
			for (FivePrimeUtrExon exon : fivePrimeUtr.getExons()) {
				ReferenceSequence referenceSeq = reference.getSubsequenceAt(exon.getChr(), exon.getStart(), exon.getEnd());
				String referenceString = (new String(referenceSeq.getBases())).toUpperCase();
				refBases = refBases + referenceString;
				if (overlaps == exon) {
					// Need to apply variant.
					String altExon = referenceString.substring(0, pos - exon.getStart()) + alt + referenceString.substring(pos - exon.getStart() + ref.length());
					altBases = altBases + altExon;
				} else {
					altBases = altBases + referenceString;
				}
			}
		} else {
			for (FivePrimeUtrExon exon : fivePrimeUtr.getExons()) {
				ReferenceSequence referenceSeq = reference.getSubsequenceAt(exon.getChr(), exon.getStart(), exon.getEnd());
				String referenceString = (new String(referenceSeq.getBases())).toUpperCase();
				refBases = refBases + reverse(referenceString);
				if (overlaps == exon) {
					// Need to apply variant.
					String altExon = referenceString.substring(0, pos - exon.getStart()) + alt + referenceString.substring(pos - exon.getStart() + ref.length());
					altBases = altBases + reverse(altExon);
				} else {
					altBases = altBases + reverse(referenceString);
				}
			}
		}
		// Find the ORFs in both versions of the 5-prime UTR
		List<Uorf> refUorfs = new ArrayList<Uorf>();
		List<Uorf> altUorfs = new ArrayList<Uorf>();
		String refBase0 = findUorfs(refBases, refBases.length() % 3, refUorfs);
		String altBase0 = findUorfs(altBases, altBases.length() % 3, altUorfs);
		String refBase1 = findUorfs(refBases, (refBases.length() + 1) % 3, refUorfs);
		String altBase1 = findUorfs(altBases, (altBases.length() + 1) % 3, altUorfs);
		String refBase2 = findUorfs(refBases, (refBases.length() + 2) % 3, refUorfs);
		String altBase2 = findUorfs(altBases, (altBases.length() + 2) % 3, altUorfs);
		String[] visualisations = new String[] {refBase0, altBase0, refBase1, altBase1, refBase2, altBase2};
		// Sort the ORF lists by consequence. The most "damaging" ORF will be first in each list
		Collections.sort(refUorfs);
		Collections.sort(altUorfs);
		// Find the most "damaging" ORF in the reference and the alternate alleles
		Uorf refUorf = null;
		Uorf altUorf = null;
		if (!refUorfs.isEmpty()) {
			refUorf = refUorfs.get(0);
		}
		if (!altUorfs.isEmpty()) {
			altUorf = altUorfs.get(0);
		}
		// This is the length change of the InDel, so we can allow for changes in distance from the gene
		int lengthChange = Math.abs(ref.length() - alt.length());
		// Work out the consequence of the change
		boolean altStats = true;
		String effect = "No change";
		if (refUorf == null) {
			if (altUorf == null) {
				return new UorfResult("", false, null, null, null, visualisations);
			} else if (altUorf.getType() == UorfType.FRAMESHIFT) {
				effect = "out-of-frame_oORF";
			} else if (altUorf.getType() == UorfType.EXTENDING) {
				effect = "CDS_elongated";
			} else {
				effect = "uORF_created";
			}
		} else if (refUorf.getType() == UorfType.FRAMESHIFT) {
			if (altUorf == null) {
				effect = "loss_out-of-frame_oORF";
				altStats = false;
			} else if (altUorf.getType() == UorfType.FRAMESHIFT) {
				if (altUorf.getStrength() > refUorf.getStrength()) {
					effect = "stronger_out-of-frame_oORF";
				} else if (altUorf.getStrength() < refUorf.getStrength()) {
					effect = "loss_weaker_out-of-frame_oORF";
					altStats = false;
				} else if (altUorf.getDistance() < refUorf.getDistance() - lengthChange) {
					effect = "closer_out-of-frame_oORF";
				} else if (altUorf.getDistance() > refUorf.getDistance() + lengthChange) {
					effect = "loss_further_out-of-frame_oORF";
					altStats = false;
				}
			} else {
				effect = "loss_out-of-frame_oORF";
				altStats = false;
			}
		} else if (refUorf.getType() == UorfType.EXTENDING) {
			if (altUorf == null) {
				effect = "loss_CDS_elongated";
				altStats = false;
			} else if (altUorf.getType() == UorfType.FRAMESHIFT) {
				effect = "out-of-frame_oORF";
			} else if (altUorf.getType() == UorfType.EXTENDING) {
				if (altUorf.getStrength() > refUorf.getStrength()) {
					effect = "stronger_CDS_elongated";
				} else if (altUorf.getStrength() < refUorf.getStrength()) {
					effect = "loss_weaker_CDS_elongated";
					altStats = false;
				} else if (altUorf.getDistance() < refUorf.getDistance() - lengthChange) {
					effect = "closer_CDS_elongated";
				} else if (altUorf.getDistance() > refUorf.getDistance() + lengthChange) {
					effect = "loss_further_CDS_elongated";
					altStats = false;
				}
			} else {
				effect = "loss_CDS_elongated";
				altStats = false;
			}
		} else {
			if (altUorf == null) {
				effect = "loss_uORF";
				altStats = false;
			} else if (altUorf.getType() == UorfType.FRAMESHIFT) {
				effect = "out-of-frame_oORF";
			} else if (altUorf.getType() == UorfType.EXTENDING) {
				effect = "CDS_elongated";
			} else {
				if (altUorf.getStrength() > refUorf.getStrength()) {
					effect = "stronger_uORF_created";
				} else if (altUorf.getStrength() < refUorf.getStrength()) {
					effect = "loss_weaker_uORF";
					altStats = false;
				} else if (altUorf.getDistance() < refUorf.getDistance() - lengthChange) {
					effect = "closer_uORF_created";
				} else if (altUorf.getDistance() > refUorf.getDistance() + lengthChange) {
					effect = "loss_further_uORF";
					altStats = false;
				}
			}
		}
		return new UorfResult(effect, !altStats, altStats ? altUorf : refUorf, refUorfs, altUorfs, visualisations);
	}

	/**
	 * Class containing the result from a uORF search.
	 */
	public static class UorfResult
	{
		private String effect;
		private boolean loss;
		private Uorf uorf;
		private List<Uorf> refUorfs, altUorfs;
		private String[] visualisations;

		public UorfResult(String effect, boolean loss, Uorf uorf, List<Uorf> refUorfs, List<Uorf> altUorfs, String[] visualisations) {
			this.effect = effect;
			this.loss = loss;
			this.uorf = uorf;
			this.refUorfs = refUorfs;
			this.altUorfs = altUorfs;
			this.visualisations = visualisations;
		}

		/**
		 * Returns a text description of the effect on uORFs of a variant.
		 *
		 * @return a String
		 */
		public String getEffect() {
			return effect;
		}

		/**
		 * Returns true if the uORF change should be considered a weakening of uORF effect, for example the variant removing a uORF.
		 *
		 * @return a boolean
		 */
		public boolean isLoss() {
			return loss;
		}

		/**
		 * Returns the uORF object that is the most relevant for the effect of the variant. If the uORF is a loss, then this uORF will be the uORF from the reference that was weakened or removed. Otherwise, this will be the uORF that was strengthened or created.
		 *
		 * @return a Uorf
		 */
		public Uorf getUorf() {
			return uorf;
		}

		/**
		 * Returns the full list of uORFs found in the reference 5-prime UTR.
		 *
		 * @return a List of Uorf objects
		 */
		public List<Uorf> getRefUorfs() {
			if (refUorfs == null) {
				return null;
			}
			return Collections.unmodifiableList(refUorfs);
		}

		/**
		 * Returns the full list of uORFs found in the alternate 5-prime UTR.
		 *
		 * @return a List of Uorf objects
		 */
		public List<Uorf> getAltUorfs() {
			if (altUorfs == null) {
				return null;
			}
			return Collections.unmodifiableList(altUorfs);
		}

		/**
		 * Returns a visualisation of the uORFs found in the reference and alternate 5-prime UTRs.
		 * This is a six-element array. Positions 0, 2, and 4 refer to the reference 5-prime UTR, and positions 1, 3, and 5 refer to the alternate 5-prime UTR. The three visualisations show the uORFs found in the three coding frames.
		 * Each String shows the bases in the 5-prime UTR, separated by spaces into codons. ORFs are shown by capitalising the bases. The number replacing the space just after the ATG codon at the beginning of an ORF is the strength of the start codon, from 1 (weak) to 3 (strong).
		 *
		 * @return an array of Strings
		 */
		public String[] getVisualisations() {
			return visualisations;
		}
	}

	private static String findUorfs(String bases, int offset, List<Uorf> uorfs) {
		StringBuilder visualisation = new StringBuilder();
		if (offset > 0) {
			visualisation.append(bases.substring(0, offset).toLowerCase() + " ");
		}
		boolean inUorf = false;
		int strength = 0;
		int distance = 0;
		int i = 0;
		for (i = offset; i < bases.length() - 2; i += 3) {
			String codon = bases.substring(i, i + 3);
			if (inUorf) {
				visualisation.append(codon + " ");
				if ("TAA".equals(codon) || "TAG".equals(codon) || "TGA".equals(codon)) {
					uorfs.add(new Uorf(distance, bases.length()  + 3 - i, strength, UorfType.NON_OVERLAPPING));
					inUorf = false;
					strength = 0;
					distance = 0;
				}
			} else {
				if ("ATG".equals(codon)) {
					inUorf = true;
					strength = 1;
					if (i >= 3) {
						if ((bases.charAt(i - 3) == 'A') || (bases.charAt(i - 3) == 'G')) {
							strength++;
						}
					}
					if (i + 3 < bases.length()) {
						if ((bases.charAt(i + 3) == 'G')) {
							strength++;
						}
					}
					distance = bases.length() - i;
					visualisation.append(codon + strength);
				} else {
					visualisation.append(codon.toLowerCase() + " ");
				}
			}
		}
		if (inUorf) {
			uorfs.add(new Uorf(distance, 0, strength, (i == bases.length() ? UorfType.EXTENDING : UorfType.FRAMESHIFT)));
		}
		visualisation.append(inUorf ? bases.substring(i) : bases.substring(i).toLowerCase());
		return visualisation.toString();
	}

	private static String reverse(String bases) {
		StringBuilder retval = new StringBuilder();
		for (int i = bases.length() - 1; i >= 0; i--) {
			char c = bases.charAt(i);
			if (c == 'A') {
				retval.append("T");
			} else if (c == 'C') {
				retval.append("G");
			} else if (c == 'G') {
				retval.append("C");
			} else if (c == 'T') {
				retval.append("A");
			} else {
				throw new RuntimeException("Invalid base sequence " + bases);
			}
		}
		return retval.toString();
	}

	private int distance;
	private int stopDistance;
	private int strength;
	private UorfType type;

	public Uorf(int distance, int stopDistance, int strength, UorfType type) {
		this.distance = distance;
		this.stopDistance = stopDistance;
		this.strength = strength;
		this.type = type;
	}

	/**
	 * Returns the distance of the start codon from the start of the gene
	 *
	 * @return an int
	 */
	public int getDistance() {
		return distance;
	}

	/**
	 * Returns the number of bases in-between the end of the ORF and the start of the gene
	 *
	 * @return an int
	 */
	public int getStopDistance() {
		return stopDistance;
	}

	/**
	 * Returns the strength of the start codon, as a number from 1 (weak) to 3 (strong)
	 *
	 * @return an int
	 */
	public int getStrength() {
		return strength;
	}

	/**
	 * Returns the strength of the start codon, as a text description
	 *
	 * @return a String
	 */
	public String getStrengthString() {
		return (strength == 1 ? "Weak" : (strength == 2 ? "Medium" : "Strong"));
	}

	/**
	 * Returns the type of the ORF.
	 *
	 * @return a UorfType
	 */
	public UorfType getType() {
		return type;
	}

	/**
	 * Returns a textual description of the ORF.
	 *
	 * @return a String
	 */
	public String toString() {
		return "(" + type + ", distance: " + distance + " to " + stopDistance + (strength == 1 ? ", Weak)" : (strength == 2 ? ", Medium)" : ", Strong)"));
	}

	/**
	 * Compare this ORF with another one, based on consequence.
	 * If this ORF is more "damaging" than the parameter, then a negative number is returned.
	 *
	 * @return an int
	 */
	public int compareTo(Uorf u) {
		if ((type == UorfType.FRAMESHIFT) && ((u.type == UorfType.EXTENDING) || (u.type == UorfType.NON_OVERLAPPING))) {
			return -1;
		} else if ((type == UorfType.EXTENDING) && (u.type == UorfType.NON_OVERLAPPING)) {
			return -1;
		} else if ((type == UorfType.EXTENDING) && (u.type == UorfType.FRAMESHIFT)) {
			return 1;
		} else if ((type == UorfType.NON_OVERLAPPING) && ((u.type == UorfType.EXTENDING) || (u.type == UorfType.FRAMESHIFT))) {
			return 1;
		} else if (strength > u.strength) {
			return -1;
		} else if (strength < u.strength) {
			return 1;
		} else if (distance < u.distance) {
			return -1;
		} else if (distance > u.distance) {
			return 1;
		} else {
			return 0;
		}
	}

	public enum UorfType
	{
		NON_OVERLAPPING, EXTENDING, FRAMESHIFT
	}
}
