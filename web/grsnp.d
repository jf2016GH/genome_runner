upstream apps {
   server 127.0.0.1:8080;
   server 127.0.0.1:8081;
}

gzip_http_version 1.0;
gzip_proxied      any;
gzip_min_length   500;
gzip_disable      "MSIE [1-6]\.";
gzip_types        text/plain text/xml text/css
                  text/javascript
                  application/javascript;

server {
   client_max_body_size 10M;
   listen 80;
   #server_name  http://162.216.114.51;

   access_log  /home/mdozmorov/logs/nginx.grsnp.log combined;
   error_log  /home/mdozmorov/logs/nginx.grsnp.log;
   
   location / {
      proxy_pass         http://127.0.0.1:8080; # port number for gr-server
      #proxy_redirect     off;
      #proxy_set_header   Host $host;
      #proxy_set_header   X-Real-IP $remote_addr;
      #proxy_set_header   X-Forwarded-For $proxy_add_x_forwarded_for;
      #proxy_set_header   X-Forwarded-Host $server_name;
   }

   location ^~ /shiny-gr {
      proxy_pass         http://127.0.0.1:3838; # port number for shiny process
      #proxy_redirect     off;
      #proxy_set_header   Host $host;
      #proxy_set_header   X-Real-IP $remote_addr;
      #proxy_set_header   X-Forwarded-For $proxy_add_x_forwarded_for;
      #proxy_set_header   X-Forwarded-Host $server_name;
   }

}
