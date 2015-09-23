Important configuration files
===
- `/etc/nginx/nginx.conf` - global nginx configuration file
- `/etc/nginx/sites-available/grsnp.d` - local nginx configuration
- `/etc/shiny-server/shiny-server.conf` - shiny server configuration

Configure nginx
===

Locale errors are caused by locale not being set on the server.

    sudo locale-gen en_US.UTF-8

Install nginx and follow [these instructions](https://www.digitalocean.com/community/tutorials/how-to-install-nginx-on-ubuntu-12-04-lts-precise-pangolin).

    sudo apt-get install nginx

There is useful information on [configuring nginx](https://www.digitalocean.com/community/tutorials/how-to-deploy-cherrypy-web-applications-behind-nginx-reverse-proxy).

Remove nginx's defaults, they lock up port 80 and intercept our calls.

    sudo rm /etc/nginx/sites-enabled/default
    sudo rm /etc/nginx/sites-available/default

Setup to run on startup:

    update-rc.d nginx defaults

Get the inet addr (i.e http://162.216.114.51/). We will use this address in the grsnp.d file.

    ifconfig

Create `grsnp.d` config file, please read it and note the comments. Edit 'access_log' and 'error_log' paths in this file. Create corresponding folders, e.g., '/home/ubuntu/logs/'. Then softlink with the sites-enabled directory.

    sudo vim /etc/nginx/sites-available/grsnp.d
    ln -s /etc/nginx/sites-available/grsnp.d /etc/nginx/sites-enabled/grsnp.d

Another configuration file: /etc/nginx/nginx.conf. You can use 

    sudo nginx -t

to check for syntax errors in the config files

    sudo service nginx start
    sudo service nginx status

Check errors

    sudo tail /var/log/nginx/error.log

Configuring R.Genomerunner
===
Next, we need to install R.genomerunner into the /srv/shiny-server folder. All folders in this /srv/shiny-server/ folder are treated as apps. Here we soft link “R.genomeruner” as “shiny-gr”. IMPORTANT the folder given to R.genomerunner in /srv/shiny-server is important. This is the name we will use for the url in results_shiny.mako.

    ln -s /home/mdozmorov/R.genomerunner/ /srv/shiny-server/shiny-gr

Open genome_runner/grsnp/frontend/templates/results_shiny.mako and edit line 98 to be:

    <iframe id="example1" style="border: none;height: 1000px; width: 100%" src="http://162.216.114.51:shiny-gr?id=${run_id}"frameborder="0"></iframe>

Critical part is to replace 'src="http://10.0.2.15:4494?id=${run_id}' by 'src="http://162.216.114.51:shiny-gr?id=${run_id}"'. 

IMPORTANT Be sure to change the file folder permission for the database. Otherwise, shiny will crash when trying to create the heatmap as it needs access to the rds file generated in the results folder.

    chown -R username:shiny-apps /home/mdozmorov/db_5.00_07-22-2015/

Also, be sure to change the `grsnp.d` config file as well so that it re-routs requests correctly to the shiny-server. i.e do: location ^~ /shiny-gr

Configure shiny-server
===

[Good overview](http://deanattali.com/2015/05/09/setup-rstudio-shiny-server-digital-ocean/#install-shiny) of installing shiny-server and user permissions.

If we need to start/stop the server.

    sudo start shiny-server
    sudo stop shiny-server
    sudo restart shiny-server

Note that you need to do reload if you want to re-read config file for shiny-server. Restart will NOT reload the config file

    sudo reload shiny-server
    sudo status shiny-server

Config file for shiny is located at 

    sudo vim /etc/shiny-server/shiny-server.conf

Read it and look for the port number located on line 5 (i.e. ‘listen 3838’)
This needs to be put into our `grsnp.d` config file created earlier.

Log file information for shiny app can be found by going to /var/log/shiny-server/. This will refresh automatically when the log file changes

    sudo tail -f shiny-gr-shiny-####-####-####.log

Other notes
===

Located the problem caused by "Cannot open Rplots.pdf".

(function (file = ifelse(onefile, "Rplots.pdf", "Rplot%03d.pdf"),  : 
  cannot open file 'Rplots.pdf')
which is the default device driver. See below

To fix:

    http://yihui.name/en/2010/12/a-special-graphics-device-in-r-the-null-device/
    Essentially, we add the code to server.R
    options(device = function(...) {
        .Call("R_GD_nullDevice", PACKAGE = "grDevices")
    })

[https://github.com/rstudio/shiny-server/issues/96](https://github.com/rstudio/shiny-server/issues/96)

Shortcuts
===

Start the server (note the port number should correspond to the `grsnp.d` settings)

    gr-server -d /home/mdozmorov/db_5.00_07-22-2015/ -r /home/mdozmorov/db_5.00_07-22-2015/ -g hg19 -p 8080

Stop all workers

    ps aux | grep "celery worker" | awk '{print $2}' | xargs kill -9

Use celery to check active/manually start workers

    celery inspect active --broker redis://localhost:7775/0
    celery worker --app=grsnp.worker_hypergeom4 -d /home/mdozmorov/db_5.00_07-22-2015/ -r /home/mdozmorov/db_5.00_07-22-2015/ --loglevel INFO -E

Troubleshooting
===
Check if something is running on port 80

    lsof -i :80 or sudo fuser -k 80/tcp 
    sudo netstat -pan | grep ":80"

Force stop running nginx

    sudo fuser -k 80/tcp

Check what listens

    sudo grep -r Listen /etc

Dealing with encoding issues

    iconv -t UTF-8 -c gf_descriptions.test > gf_descriptions.test.utf8

413 Request Entity Too Large - [check this link](http://cnedelcu.blogspot.com/2013/09/nginx-error-413-request-entity-too-large.html) and add 'client_max_body_size 10M;' to the nginx config file

ERROR parallel: Warning: No more processes: Decreasing number of running jobs to 1. Raising ulimit -u or /etc/security/limits.conf may help.

    ulimit -u # Check the allowed limit of concurrent processes
    ps -u <username> # Count the number of running processes

[Linux system monitoring tools](http://www.cyberciti.biz/tips/top-linux-monitoring-tools.html)
