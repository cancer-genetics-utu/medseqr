byte.convert <- function
		(x, to=c("b", "k", "m", "g", "t"), digits=2)
{
	ptn <- "(\\d*(.\\d+)*)(.*)"
	num  <- as.numeric(sub(ptn, "\\1", x))
	unit <- tolower(sub(ptn, "\\3", x))
	unit[unit==""] <- "b"
	mult <- c("b"=1, "k"=1024, "m"=1024^2, "g"=1024^3, "t"=1024^4)
	if (is.null(to) && all(unit=="b")) {
		to <- names(mult)[ceiling(nchar(x)/3)]
	} else {
		to <- match.arg(to)
	}
	if (all(unit=="b")) {
		paste0(round(num / unname(mult[to]), digits), " ", ifelse(to=="k", to, toupper(to)), ifelse(to=="b", "", "B"))
	} else if (to!="b") {
		byt <- num * unname(mult[unit])
		paste0(round(byt / unname(mult[to]), digits), " ", ifelse(to=="k", to, toupper(to)), "B")
	} else {
		num * unname(mult[unit])
	}
}

