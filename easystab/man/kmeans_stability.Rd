\name{kmeans_stability}

\alias{kmeans_stability}

\title{Adapts a list of clusterings from \code{kmeans} to one usable by
  \code{perturbationStability}.}

\description{Given a list of clusterings, each from \code{kmeans},
  returns a corresponding list of clusterings suitable for input to
  \code{perturbationStability}.}

\usage{ kmeans_stability(x, kmeans_output_list)}

\arguments{
  \item{x}{Matrix or data frame object containing the clustered data.}
  
  \item{kmeans_output_list}{A list of clusterings, each being the output of the
    kmeans function.}
}

\value{ A list of clusterings that can be used as imput to the
  \code{perturbationStability} function.}

\seealso{
  \code{\link{easystab}},
  \code{\link{perturbationStability}},
  \code{\link{hclust_stability}}
}

