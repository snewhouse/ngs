#!/bi/bash

VERSION=${1} 
TAG=${1}

DATE=`date +"%Y%m%d%H%M%S`

# add all and update dev2
git checkout dev2
git pull
git add --all ./*
git commit -am "update ${DATE}"
git push
git pull

# upate
git checkout dev2
git pull
git checkout master
git pull

# merge master with dev2
git checkout master
git merge --no-ff dev2 -m "latest updates ${DATE}"
git push
git pull

# creating release branch
git checkout master
git pull

git checkout -b ngseasy-${VERSION} master

git add --all ./*

git commit -am "ngseasy-${VERSION} ${DATE}"

git tag -a ${TAG}

git push origin ngseasy-${VERSION}

git checkout master
git branch -D ngseasy-${VERSION}

git checkout dev2
git pull

#http://nvie.com/posts/a-successful-git-branching-model/