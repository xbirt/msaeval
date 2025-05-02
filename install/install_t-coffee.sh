#!/bin/bash

#sudo apt update && sudo apt install -y --no-install-recommends t-coffee

wget https://s3.eu-central-1.amazonaws.com/tcoffee-packages/Stable/Latest/T-COFFEE_installer_Version_13.46.0.919e8c6b_linux_x64.bin
chmod 755 T-COFFEE_installer_Version_13.46.0.919e8c6b_linux_x64.bin
sudo ./T-COFFEE_installer_Version_13.46.0.919e8c6b_linux_x64.bin --mode unattended --user_email student@usv.ro --prefix /usr/local/t_coffee
sudo ln -s /usr/local/t_coffee/bin/t_coffee /usr/bin
