# loadflexBatch

This repository contains a demo script and data for running the `loadflex` 
package in batch mode, i.e., for multiple sites, constituents, and models.

While we are committed to sharing our scripts openly, and you are welcome to
try the scripts yourself, the contents of this repository were created for a
specific application, and we regret that we do not have the resources to update
these scripts or provide support for their use by others.

You will likely have the most success with these `loadflexBatch` scripts if you use
`loadflex` version 1.2.0, which can be installed with

```r
devtools::install_github('USGS-R/loadflex@v1.2.0')
```

For a description of this repository's contents, and instructions on how to use 
the scripts contained here, see 
[batch_tutorial.html](https://usgs-r.github.io/loadflexBatch/batch_tutorial.html).

Updated 7/12/2019

# Note on Branch Names
On 2022-02-16, the `master` branch of this repository was renamed from `master` to `main`, as part of a coordinated change across the USGS-R and USGS-VIZLAB organizations.  This is part of a broader effort across the git community to use more inclusive language. For instance, [git](https://sfconservancy.org/news/2020/jun/23/gitbranchname/), [GitHub](https://github.com/github/renaming), and [GitLab](https://about.gitlab.com/blog/2021/03/10/new-git-default-branch-name/) have all changed or are in the process of changing their default branch name. If you have forks or local clones dating to before 2022-02-16, you will probably want to update them as well.

### Changing default branches on forks and local clones
First, update your local git configuration so that new repositories created locally will have the correct default branch: `git config --global init.defaultBranch main`.

Now, for any forks, do the following:
1. Go to \<your fork\> -> Settings -> Branches and edit the default branch from `master` to `main`.
1. Update the settings for your local clone of this fork to match this change.

```r
git branch -m master main
git fetch origin
git branch -u origin/main main
git remote set-head origin -a
```
Finally, search within your repository for `master` so that you can change references (e.g. URLs pointing to specific commit points of files) to point to `main` instead.
