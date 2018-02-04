def Consensus(infile, outfile, W, Tmm):
# JBB version:
# Look for 2 repeats in the sequence with max length as 600

	import sys	#JHL, 20160303
	import re	#JHL, 20160303
	import math
	from scipy.stats import mode
	from itertools import izip

	cdef int i
	cdef int counter_PoorQuality
	cdef int counter_NoRepeats
	cdef int counter_AbnormalRepeatLength
	cdef int counter_LowIdentity
	cdef int counter_TotalReads
	cdef int counter_ConsensusSequences

	cdef int StartPosition
	cdef int Identity

	counter_PoorQuality = 0
	counter_NoRepeats = 0
	counter_AbnormalRepeatLength = 0
	counter_LowIdentity = 0
	counter_TotalReads = 0
	counter_ConsensusSequences = 0
	RepeatLengths = [0]*600 #JBB, 20170808. Save repeat length histogram

	#W = 12 # window size
	#Tmm = 2 # Threshold for mismatches at either-end extension

	while True:

		SequenceID = infile.readline() # 1st line of fastq
		if not SequenceID: # at the end of the fastq file
			break

		counter_TotalReads += 1
		Sequence = infile.readline() # 2nd line of fastq
		Sequence = Sequence.rstrip("\n")
		EmptyLine = infile.readline() # 3rd line of fastq: +
		QualityScores = infile.readline() # 4th line of fastq
		QualityScores = QualityScores.rstrip("\n")

		#Remove the first base which can be error prone
		Sequence = Sequence[1:]
		ReadLength = len(Sequence)
		QualityScores = QualityScores[1:]

		if Sequence.count("N") > 0.1 * ReadLength: #JBB, 20170808: "90" -> "0.1 * ReadLength"
			counter_PoorQuality += 1

		else:
			L = 0 # the length of the 2nd string in Substrings
			i = 0
			while i <= ReadLength - (W + 1):
				Substrings = Sequence.split(Sequence[i:i+W]) #JBB, 20170808, win size is flexible, depending on the data quality
				numStrings = len(Substrings)
				if numStrings >= 3:
					L = len(Substrings[1])
					break
				i += 1

			#Remove/count reads with no detectable repeats
			if L == 0:
				counter_NoRepeats += 1

			elif L < W/2 or L > ReadLength - 2*W:
				counter_AbnormalRepeatLength += 1

			else:
				# Extend towards the both ends to find the repeats.
				u1 = i # upstream index of R1
				d1 = i + W - 1 # downstream index of R1
				R1 = Sequence[u1:d1]

				u2 = i + W + L
				d2 = i + W + L + W - 1
				R2 = Sequence[u2:d2]

				umm = 0 # mismatch counter at upstream
				dmm = 0

				while umm < Tmm and u1 > 0:
					R1 = Sequence[u1-1] + R1
					R2 = Sequence[u2-1] + R2
					if Sequence[u1-1] != Sequence[u2-1]:
						umm += 1
					u1 -= 1
					u2 -= 1
				while dmm < Tmm and d2 < ReadLength and d1 < u2:
					R1 = R1 + Sequence[d1]
					R2 = R2 + Sequence[d2]
					if Sequence[d1] != Sequence[d2]:
						dmm += 1
					d1 += 1
					d2 += 1

				tt = Tmm
				while tt:
					if R1[0] != R2[0]:
						R1 = R1[1:]
						R2 = R2[1:]
						u1 += 1
						u2 += 1
					if R1[-1] != R2[-1]:
						R1 = R1[:-1]
						R2 = R2[:-1]
						d1 -= 1
						d2 -= 1
					tt -= 1

				RepeatLength = len(R1)
				RepeatLengths[RepeatLength] += 1

				if RepeatLength < 2*W:
					counter_LowIdentity += 1
				else:
					counter_ConsensusSequences += 1
					ConsensusSequence = ""
					ReCalculatedQualityScores = ""

					for i in range(RepeatLength):
						Q1 = ord(QualityScores[u1]) # ord() returns ASCII code.
						Q2 = ord(QualityScores[u2])
						if R1[i] == R2[i] and R1[i] != "N":
							ConsensusSequence += R1[i]
							ReCalculatedQualityScores += chr(int(round((Q1 + Q2)/2.0)))
						else:
							if Q1 >= Q2:
								ConsensusSequence += R1[i]
								ReCalculatedQualityScores += QualityScores[u1]
							else:
								ConsensusSequence += R2[i]
								ReCalculatedQualityScores += QualityScores[u2]
						u1 += 1
						u2 += 1

					outfile.write(SequenceID)
					outfile.write(ConsensusSequence + "\n")
					outfile.write(EmptyLine)
					outfile.write(ReCalculatedQualityScores + "\n")

	return counter_PoorQuality, counter_NoRepeats, counter_AbnormalRepeatLength, counter_LowIdentity, counter_ConsensusSequences, counter_TotalReads, RepeatLengths
