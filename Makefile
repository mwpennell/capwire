all: install

document: roxygen staticdocs

roxygen:
	@mkdir -p man
	Rscript -e "library(methods); devtools::document()"

staticdocs:
	@mkdir -p inst/staticdocs
	Rscript -e "library(methods); staticdocs::build_site()"

publish_pages:
	cd inst && ./update-gh-pages.sh

install:
	R CMD INSTALL --no-test-load .

build:
	R CMD build .

check: build
	R CMD check --no-manual `ls -1tr capwire*gz | tail -n1`
	@rm -f `ls -1tr capwire*gz | tail -n1`
	@rm -rf capwire.Rcheck

# No real targets!
.PHONY: all test document install
