\name{fastqc.summary}
\alias{fastqc.summary}
\title{Generate summary statistics from FastQC output.}
\usage{
  fastqc.summary(files, outfolder = NULL, out.id = NULL, top = 5,
                 pdf.width = 12, text.size = 9, title = NULL)
}
\arguments{
  \item{files}{\code{character}. Vector of FastQC output summary file paths (summary.txt).}

  \item{outfolder}{\code{character}. Optional output folder.}

  \item{out.id}{\code{character}. Optional output identifier. Is appended to the file name.}

  \item{top}{\code{integer}. Number of files for which the summary statistics will be returned, sorted by quality (worst to best).}

  \item{pdf.width}{\code{integer}. Width of PDF in inches.}

  \item{text.size}{\code{numeric}. Size of text elements in plots.}

  \item{title}{\code{character}. Optional amendment plot title.}
}
\value{
  A list of summary statistics with three slots: 'per_sample', 'per_variable' and 'all_param_per_sample'
}
\description{
  Generate summary statistics from FastQC output.
}
\details{
  Given a vector of file paths to fastqc summary.txt-files, returns statistics for the files.
  Encodes fastqc scale PASS-WARN-FAIL to numeric from 0 to 2.
  Optionally, three summary plots are produced.
}
\seealso{
  \code{\link{ggplot}}
}
