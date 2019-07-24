.PHONY: all
all:
	echo 'type "make doc" for updating documentation'
	echo 'then, "git commit" and "git push"'

.PHONY: doc
doc:
	rm -rf docs
	mkdir docs
	cd documentation; make clean html; cd ..
	cp -a documentation/_build/html/* docs
	touch docs/.nojekyll
	echo 'now:'
	echo 'git commit -m "updated documentation"'
	echo 'git push'
