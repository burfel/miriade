---
layout: default
title: Index
order: 1
---

# Gitlab pages
Welcome to an example page! Below you will find a detailed instruction on how to set-up your own Jekyll website with Gitlab-CI and our LCSB template.

Sources for this page are [available in Gitlab](https://git-r3lab.uni.lu/core-services/pages-jekyll-lcsb-template).
If you are interested, have found any issues with the layout or have a helpful suggestion, you can also navigate to the [repository for theme](https://git-r3lab.uni.lu/core-services/jekyll-theme-lcsb-default).

# Setting up your web page

## The whole process to deploy your website
0. Make sure that you have access to our [LCSB's Gitlab - https://git-r3lab.uni.lu/](https://git-r3lab.uni.lu/). If you ever cloned a repository or set-up a new one there, then you should be fine.  If you don't have it, either contact the sysadmins (`lcsb-sysadmins (at) uni.lu`) or open a ticket at [https://service.uni.lu](https://service.uni.lu). 
1. [Create an empty repository in Gitlab](https://git-r3lab.uni.lu/projects/new). Please bear in mind, that the both names of the namespace and the project influence 
    final address of the page - it will internally follow the `https://<namespace>.pages.uni.lu/<project_name>` convention. (Note, that in the very last step, SIU can set a new alias/URL for your website)
    ![image](assets/screenshots/new_project.png)
2. On your computer, clone this very repository using: `git clone ssh://git@git-r3lab-server.uni.lu:8022/core-services/pages-jekyll-lcsb-template.git`.
3. Once cloned, navigate into cloned repository (`cd pages-jekyll-lcsb-template`) and remove the _remote_ - `git remote rm origin` (so that you update your repository, and not this very page).
4. Set a new remote (so that you push to your new repository) - `git remote add origin ssh://git@git-r3lab-server.uni.lu:8022/firstname.surname/your_projects_name.git`. You can find the correct remote address in Gitlab just after creating new repository, as in the following image:
 ![image](assets/screenshots/remote.png)
5. Modify site's settings (`_config.yml`) to match your needs. Refer to the next section for help.
6. Modify the index page. Modify (or delete) help and about pages. Add your own content. 
7. Add your changes (`git add .`), commit (`git commit -m "Initial commit"`) and push (`git push --set-upstream origin master`) to the repository.
8. Your page is published! Go to `https://<namespace>.pages.uni.lu/name-of-repository` in your favourite browser, or the URL you specified in the SIU ticket.
9. In gitlab, go to **Settings** (under left-hand menu) > **General** > **Advanced** (hit `Expand` button) > **Remove fork relationship** (red button), then follow the instructions from the pop-up.
10. If you want to have your page publicly available - contact us  (`lcsb-sysadmins (at) uni.lu`), we will make a ticket to SIU.

## What should you change in settings file?
We used to require a change in `url` and `baseurl` - but not anymore :)

However, you still might want to change:

 * `title` field
 * `e-mail` field
 * `description` field
 * `date` field
 * `banner` field - if you want to have your own banner (the text next to _uni.lu_ logo), please contact us.
 
## Testing the web page locally
You can test your website locally (on your machine). 

* First, make sure that you have Ruby installed. If not - please [install it](https://www.ruby-lang.org/en/downloads/).
* Then, from terminal install _bundler_ - `gem install bundler`. 
* Navigate into a directory with your website
* Initialize the site with: `bundle install`.
* Finally, run the site: `bundle exec jekyll serve`.


## Please don't change/remove Gemfile and .gitlab-ci.yml files
They are mandatory in order for the website to work. First one contains the website dependencies, the second is responsible for building the website.

### If you want to have the banner link to different target (e.g. uni.lu's index site)
Uncomment `banner_link: https://wwwen.uni.lu/` in the `_config.yml`, and change it to your target URI.

## Common problems
### *The website is not updated after commiting to the repository!*
Did you push the commit? If yes, then you probably changed/deleted `.gitlab-ci.yml` file.

### *The website looks broken! There are no images, no colors etc.*
You probably didn't configure `baseurl` parameter in the settings or configured it wrongly. Please take a look on `_settings.yml` file.

### *The links in the menu are not working (they point to "404: Not found").*
You probably didn't add `permalink` attribute. Or the post has `published: false` or `draft: true` set. Please take a look on the post file.

### *Something goes wrong with Gitlab-CI*
It never happened before, please notify us immediately.

### Other problems?
Please send us an email! (`lcsb-sysadmins (at) uni.lu`)
