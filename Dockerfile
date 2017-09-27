FROM r-base

MAINTAINER Tony Cannistra version:0.1

# mount repo (must be local and run with -v $(PWD):/projdir)
VOLUME /projdir
WORKDIR /projdir

# install dependencies
RUN apt-get update && apt-get install -y libgeos-dev libgdal-dev default-jdk libssl-dev git
COPY INSTALL.r INSTALL.r
RUN Rscript INSTALL.r

# download + install maxent
ADD http://biodiversityinformatics.amnh.org/open_source/maxent/maxent.php?op=download /usr/local/lib/R/site-library/dismo/java/
RUN cd /usr/local/lib/R/site-library/dismo/java/ && unzip maxent.zip 

CMD "/usr/bin/bash"