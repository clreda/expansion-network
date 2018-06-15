crms <- list(nodes=c("S1", "S2", "G1", "G2"),
		isTF=list(
			S1=TRUE,
			S2=TRUE,
			G1=FALSE,
			G2=TRUE
			),
			crms=list(
				S1=NULL,
				S2=NULL,
				G1=list("M1", "M5"),
				G2=list("M2", "M3", "M4")
			),
			tfs=list(
				S1=NULL,
				S2=NULL,
				G1=list(M1=list("S1", "S2"), M5=list("G2")),
				G2=list(M2=list("S1"), M3=list("S2"), M4=list("S2"))
			)
		)
