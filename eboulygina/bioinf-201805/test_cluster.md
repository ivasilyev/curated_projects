> **DISCLAIMER**
> 
> ALL SCRIPTS HAVE BEEN TESTED ON 3 VMWARE WORKSTATION VMS SUPPLIED WITH UBUNTU SERVER 16.04
> 
> THERE ARE NO GURANTEES THAT EVERYTHING SHALL WORK EVEN IN THE SAME CASE

# Preparation
## Used software
- [VMware Workstation 14](https://www.vmware.com/products/workstation-pro/workstation-pro-evaluation.html)
- [Ubuntu Server 16.04](https://www.ubuntu.com/download/server)
- [YUMI – Multiboot USB Creator](https://www.pendrivelinux.com/yumi-multiboot-usb-creator/)
- [KiTTY, a lightweight telnet and SSH client for Windows](http://www.9bis.net/kitty/?page=Download)
- [Xming X Server for Windows](https://sourceforge.net/projects/xming/files/latest/download)
- [Fonts for Xming X Server for Windows](https://sourceforge.net/projects/xming/files/Xming-fonts/7.7.0.10/Xming-fonts-7-7-0-10-setup.exe/download)

To make SSH X11 forwarding work, in KiTTY go to *Connection - SSH - X11* and check *Enable X11 forwarding*. Then launch the XLaunch server.

## Cluster creation
### HOST OPERATIONS
#### Update cache on the host node
```
sudo apt-get -y update; sudo apt-get -y upgrade; sudo apt-get -y autoremove
```
#### Install packages
```
sudo apt-get -y install openssh-client openssh-server python sshpass apt-transport-https ca-certificates curl software-properties-common git python-pip python3-pip
```
##### (Optional) Install desktop environment
```
sudo apt-get -y install xubuntu-desktop gksu leafpad synaptic gnome-software terminator gedit geany
```
##### (Optional) For VMware VM with DE run:
```
sudo apt-get -y install open-vm-tools-desktop
```
#### Reboot
```
sudo shutdown -r now
```
#### View network interfaces
```
ip addr show
```
#### Scan LAN if required
```
sudo apt-get -y install arp-scan
sudo arp-scan --localnet
```
#### Add nodes to system hosts
```
sudo cp /etc/hosts /etc/hosts.bak
sudo nano /etc/hosts
```
```
127.0.0.1 localhost localhost.localdomain
# Not changing the default hostname to replicate it
127.0.1.1 ubuntu

# The following lines are desirable for IPv6 capable hosts
::1     localhost ip6-localhost ip6-loopback
ff02::1 ip6-allnodes
ff02::2 ip6-allrouters

# master
10.0.0.100 node0 node0.cluster
# minions
10.0.0.101 node1 node1.cluster
10.0.0.102 node2 node2.cluster
```
#### Change host name if required
```
sudo nano /etc/hostname
```
#### Configure SSH daemon
##### Host:
```
sudo cp /etc/ssh/sshd_config /etc/ssh/sshd_config.bak
sudo nano /etc/ssh/sshd_config
```
```
Port 2121
Protocol 2
ListenAddress 127.0.0.1
ListenAddress 10.0.0.100
# Add external ListenAddress here
SyslogFacility AUTHPRIV
LogLevel INFO
LoginGraceTime 90
PermitRootLogin no
AllowUsers user
MaxAuthTries 3
X11Forwarding yes
X11DisplayOffset 10
X11UseLocalhost yes
PrintMotd no
PrintLastLog yes
PubkeyAuthentication yes
AuthorizedKeysFile	.ssh/authorized_keys .ssh/authorized_keys2
PasswordAuthentication yes
ChallengeResponseAuthentication no
GSSAPIAuthentication yes
GSSAPICleanupCredentials yes
UsePAM yes
AcceptEnv LANG LC_*
# Look for Subsystem correct path in /etc/ssh/sshd_config.bak
Subsystem sftp /usr/lib/openssh/sftp-server
```
##### Workers (on host):
```
sudo nano /etc/ssh/sshd_config_minions
```
```
Port 22
Protocol 2
# Add internal ListenAddress here to avoid access from external network
SyslogFacility AUTHPRIV
LogLevel FATAL
LoginGraceTime 90
PermitRootLogin no
AllowUsers user
MaxAuthTries 3
X11Forwarding no
PrintMotd no
PrintLastLog yes
PubkeyAuthentication yes
AuthorizedKeysFile	.ssh/authorized_keys .ssh/authorized_keys2
PasswordAuthentication yes
ChallengeResponseAuthentication no
GSSAPIAuthentication yes
GSSAPICleanupCredentials yes
UsePAM yes
AcceptEnv LANG LC_*
# Look for Subsystem correct path in /etc/ssh/sshd_config.bak
Subsystem sftp /usr/lib/openssh/sftp-server
```

#### Apply changes
```
sudo chmod a-w /etc/ssh/sshd_config.bak /etc/hosts.bak
sudo service ssh restart
```
### NFS SERVER INITIALIZATION
```
ip addr show | grep "inet "
```
> ```
> 10.0.0.100/24 ...
> ```

```
sudo apt-get -y install nfs-kernel-server nfs-common
sudo mkdir /data
sudo chmod -R 777 /data
sudo nano /etc/exports
```
```
/data *(sec=sys,rw,async,no_root_squash)
```
```
sudo /etc/init.d/nfs-kernel-server restart
```
### WORKERS OPERATIONS
#### Now logout from the master and login to _each_ node and perform the manual configuration:
```
sudo apt-get -y update; sudo apt-get -y upgrade; sudo apt-get -y autoremove
sudo apt-get -y install openssh-client openssh-server python sshpass apt-transport-https ca-certificates curl software-properties-common git python-pip python3-pip
sudo cp /etc/ssh/sshd_config /etc/ssh/sshd_config.bak
sudo cp /etc/hosts /etc/hosts.bak
sudo chmod a-w /etc/ssh/sshd_config.bak /etc/hosts.bak
```
#### Use only host ip address during manual configuration files replication
```
export GET_CFG="scp -P 2121 user@10.0.0.100"
sudo ${GET_CFG}:/etc/ssh/sshd_config_minions /etc/ssh/sshd_config
sudo ${GET_CFG}:/etc/hosts /etc/hosts
sudo service ssh restart
```
### HOST OPERATIONS
#### Install Ansible
```
sudo apt-add-repository universe
sudo apt-add-repository ppa:ansible/ansible
sudo apt-get -y update
sudo apt-get -y install ansible
ansible --version | head -n 1
```
#### Create Ansible variables
```
# Create hosts variable
echo "export AWB_HOSTS=${HOME}/.ansible/cfg/hosts" | tee -a ~/.bashrc
export AWB_HOSTS=${HOME}/.ansible/cfg/hosts
# Create scripts directory variable
echo "export AWB_DIR=${HOME}/.ansible/scripts" | tee -a ~/.bashrc
export AWB_DIR=${HOME}/.ansible/scripts
# Create remote user name variable
echo "export AWB_UN=$(whoami)" | tee -a ~/.bashrc
export AWB_UN=$(whoami)
```
#### Start using Ansible
##### Create inventory
```
mkdir -p ~/.ansible/cfg
mkdir ${AWB_DIR}
nano ${AWB_HOSTS}
```
```
[all_hosts]
node0 ip=10.0.0.100 ansible_connection=ssh ansible_user=user ansible_port=2121
node1 ip=10.0.0.101
node2 ip=10.0.0.102

[kube-master]
node0

[etcd]
node0

[kube-node]
node1
node2

[k8s-cluster:children]
kube-node
kube-master
```
##### Create main scripts to update SSH keys
```
nano ${AWB_DIR}/setup_ssh_keys_1.sh
```
```
#!/usr/bin/env bash
# setup_ssh_keys_1.sh
echo "Disable login message"
touch ~/.hushlogin
echo "Create backup"
mkdir -p ~/ssh_backups
cat ~/.ssh/authorized_keys >> ~/ssh_backups/authorized_keys
cat ~/.ssh/known_hosts >> ~/ssh_backups/known_hosts
echo "Purge SSH directory"
rm -rf ~/.ssh
mkdir ~/.ssh
echo "Generate SSH key"
ssh-keygen -t rsa -b 4096 -f ~/.ssh/id_rsa -C $(hostname)"_key" -q -N ""
eval `ssh-agent -s`
ssh-add ~/.ssh/id_rsa
echo "Launch SSH key exchange"
printf "Host node0\n  HostName node0\n  Port 2121\n  User user\n\n"| tee -a ~/.ssh/config
```
```
nano ${AWB_DIR}/setup_ssh_keys_2.sh
```
```
#!/usr/bin/env bash
# setup_ssh_keys_2.sh
# All nodes list, host is latest
NODES=( node2 node1 node0 )
for NODE in "${NODES[@]}"
    do
        echo "Update 'authorized_keys'"
        # This line is unsafe
        sshpass -p "qwerty" ssh-copy-id -i ~/.ssh/id_rsa.pub -o 'StrictHostKeyChecking=no' user@${NODE}
        # For non-Debian OS do manually: ssh-copy-id -i ~/.ssh/id_rsa.pub user@<node>
        echo "Update 'known_hosts'"
        # ssh-keygen -R ${NODE}
        echo "#" ${NODE} >> ~/.ssh/known_hosts
        echo $(ssh-keyscan -H ${NODE}) >> ~/.ssh/known_hosts
        # echo "Send public key"
        # cat ~/.ssh/id_rsa.pub | ssh user@${NODE} "mkdir -p ~/.ssh && cat >> ~/.ssh/authorized_keys"
        echo "Processed node:" ${NODE}
    done
echo "Export complete"
chmod 700 ~/.ssh
chmod 400 ~/.ssh/id_rsa
chmod 640 ~/.ssh/authorized_keys ~/.ssh/known_hosts
chmod 600 ~/.ssh/config
```
##### Create `setup_ssh_keys.sh` automation scripts
```
nano ${AWB_DIR}/purge_ssh_environment.yml
```
```
# ansible-playbook -i ${AWB_HOSTS} --user ${AWB_UN} --ask-pass ${AWB_DIR}/purge_ssh_environment.yml
- hosts: k8s-cluster:children
  tasks:
    - name: Remove files
      file:
        path: "{{ item }}"
        state: absent
      with_items:
        - ~/.ssh/authorized_keys
        - ~/.ssh/known_hosts
        - ~/authorized_keys
        - ~/known_hosts
```
```
nano ${AWB_DIR}/setup_ssh_keys_1.yml
```
```
# ansible-playbook -i ${AWB_HOSTS} --user ${AWB_UN} --ask-pass ${AWB_DIR}/setup_ssh_keys_1.yml
- hosts: k8s-cluster:children
  tasks:
    - name: Replicate script 'setup_ssh_keys.sh'
      copy:
          src: setup_ssh_keys_1.sh
          dest: ~/setup_ssh_keys.sh
    - name: Execute script 'setup_ssh_keys.sh'
      command: bash '~/setup_ssh_keys.sh'
    - name: Remove script 'setup_ssh_keys.sh'
      file:
        state: absent
        path: ~/setup_ssh_keys.sh
```
```
nano ${AWB_DIR}/setup_ssh_keys_2.yml
```
```
# ansible-playbook -i ${AWB_HOSTS} --user ${AWB_UN} --ask-pass ${AWB_DIR}/setup_ssh_keys_2.yml
- hosts: k8s-cluster:children
  tasks:
    - name: Replicate script 'setup_ssh_keys.sh'
      copy:
          src: setup_ssh_keys_2.sh
          dest: ~/setup_ssh_keys.sh
    - name: Execute script 'setup_ssh_keys.sh'
      command: bash '~/setup_ssh_keys.sh'
    - name: Remove script 'setup_ssh_keys.sh'
      file:
        state: absent
        path: ~/setup_ssh_keys.sh
```
##### Create reboot automation script
```
nano ${AWB_DIR}/ansible-reboot-and-wait.yml
```
```
# ansible-playbook -i ${AWB_HOSTS} --user ${AWB_UN} --ask-become-pass ${AWB_DIR}/ansible-reboot-and-wait.yml
- hosts: k8s-cluster:children
  become: true
  become_method: sudo
  tasks:
      # Reboot a host and wait for it to return
    - name: Send the reboot command
      shell: shutdown -r +1
      # 60 seconds before reboot
    - name: This pause is mandatory, otherwise the existing control connection gets reused!
      pause: seconds=300
      # This works with the existing ansible hosts inventory and so any custom ansible_ssh_hosts definitions are being used
    - name: Now we will run a local 'ansible -m ping' on this host until it returns.
      local_action: shell ansible -u {{ ansible_user_id }} -m ping {{ inventory_hostname }}
      register: result
      until: result.rc == 0
      retries: 30
      delay: 10
    - name: And finally, execute 'uptime' when the host is back.
      shell: uptime
    - name: Disable swap
      shell: swapoff -a
```
##### Create `/etc/hosts` replication script
```
nano ${AWB_DIR}/replicate_hosts.yml
```
```
# ansible-playbook -i ${AWB_HOSTS} --user ${AWB_UN} --ask-become-pass ${AWB_DIR}/replicate_hosts.yml
- hosts: k8s-cluster:children
  become: true
  become_method: sudo
  tasks:
    - name: Replicate hosts
      copy:
          src: /etc/hosts
          dest: /etc/hosts
    - name: Update hostname
      shell: echo {{ inventory_hostname }} > /etc/hostname
```
##### Launch SSH keys discovery
```
# Replicate hostnames, SSH and root passwords are required
ANSIBLE_HOST_KEY_CHECKING=False ansible-playbook -i ${AWB_HOSTS} --user ${AWB_UN} --ask-pass --ask-become-pass ${AWB_DIR}/replicate_hosts.yml
# Update SSH keys, only SSH password is required
ANSIBLE_HOST_KEY_CHECKING=False ansible-playbook -i ${AWB_HOSTS} --user ${AWB_UN} --ask-pass ${AWB_DIR}/setup_ssh_keys_1.yml
# Reboot cluster, SSH and root passwords are required
ANSIBLE_HOST_KEY_CHECKING=False ansible-playbook -i ${AWB_HOSTS} --user ${AWB_UN} --ask-pass --ask-become-pass ${AWB_DIR}/ansible-reboot-and-wait.yml
# Purge SSH environment
ANSIBLE_HOST_KEY_CHECKING=False ansible-playbook -i ${AWB_HOSTS} --user ${AWB_UN} --ask-pass ${AWB_DIR}/purge_ssh_environment.yml
# Validate SSH keys, only SSH password is required
ANSIBLE_HOST_KEY_CHECKING=False ansible-playbook -i ${AWB_HOSTS} --user ${AWB_UN} --ask-pass ${AWB_DIR}/setup_ssh_keys_2.yml
# Clean host
rm -f ${AWB_DIR}/setup_ssh_keys_1.sh ${AWB_DIR}/setup_ssh_keys_2.sh
```
#### Mount NAS devices (node the blank lines flanking the command)
```
nano ${AWB_DIR}/fstab
```
```

10.0.0.100:/data /data nfs rw,sec=sys,vers=4,addr=10.0.0.100,clientaddr=10.$ 0 0

```
```
nano ${AWB_DIR}/append_fstab.yml
```
```
# ansible-playbook -i ${AWB_HOSTS} --user ${AWB_UN} --ask-become-pass ${AWB_DIR}/append_fstab.yml
- hosts: kube-node
  become: true
  become_method: sudo
  tasks:
    - name: Copy fstab
      copy:
          src: fstab
          dest: ~/fstab
    - name: Update fstab
      shell: cat ~/fstab | tee -a /etc/fstab
    - name: Cleanup
      file:
        state: absent
        path: ~/fstab
```
#### Create swap disabling script to run k8s
```
nano ${AWB_DIR}/disable_swap.yml
```
```
# ansible-playbook -i ${AWB_HOSTS} --user ${AWB_UN} --ask-become-pass ${AWB_DIR}/disable_swap.yml
- hosts: k8s-cluster:children
  become: true
  become_method: sudo
  tasks:
    - name: Disable swap
      shell: swapoff -a; sleep 10s
    - name: Erase swap file
      shell: rm -rf /swapfile
    - name: Backup '/etc/fstab'
      shell: cp /etc/fstab /etc/fstab_with_swap.bak
    - name: Unmount swap
      shell: sed -i '/swap/s/^/#/' /etc/fstab
```
```
ansible-playbook -i ${AWB_HOSTS} --user ${AWB_UN} --ask-become-pass ${AWB_DIR}/append_fstab.yml
ansible-playbook -i ${AWB_HOSTS} --user ${AWB_UN} --ask-become-pass ${AWB_DIR}/disable_swap.yml
ansible-playbook -i ${AWB_HOSTS} --user ${AWB_UN} --ask-become-pass ${AWB_DIR}/ansible-reboot-and-wait.yml
```
### Other Ansible charts
#### Grant permissions for Docker (if already installed)
```
nano ${AWB_DIR}/usermod_docker.yml
```
```
# ansible-playbook -i ${AWB_HOSTS} --user ${AWB_UN} --ask-become-pass ${AWB_DIR}/usermod_docker.yml
- hosts: k8s-cluster:children
  become: true
  become_method: sudo
  tasks:
    - name: Grant privileges
      shell: usermod -aG docker user
```
#### Synchronize bash aliases
```
nano ${AWB_DIR}/replicate_bash_aliases.yml
```
```
# ansible-playbook -i ${AWB_HOSTS} --user ${AWB_UN} ${AWB_DIR}/replicate_bash_aliases.yml
- hosts: kube-node
  tasks:
    - name: Copy bash_aliases
      copy:
          src: ~/.bash_aliases
          dest: ~/.bash_aliases
```
#### Install software packages
```
nano ${AWB_DIR}/install_soft.yml
```
```
# ansible-playbook -i ${AWB_HOSTS} --user ${AWB_UN} --ask-become-pass ${AWB_DIR}/install_soft.yml
- hosts: kube-node
  become: true
  become_method: sudo
  tasks:
    - name: Manage apt
      apt:
        upgrade: yes
        update_cache: yes
    - name: Install list of packages
      apt: name={{ item }} state=present
      with_items:
       - openssh-client
       - openssh-server
       - python
       - sshpass
       - apt-transport-https
       - ca-certificates
       - curl
       - software-properties-common
       - git
       - python-pip
       - python3-pip
    - name: Clean apt
      apt:
        autoremove : yes
        autoclean: yes

```
#### Create cdrom repository disabling script
```
nano ${AWB_DIR}/disable_cdrom_repo.yml
```
```
# ansible-playbook -i ${AWB_HOSTS} --user ${AWB_UN} --ask-become-pass ${AWB_DIR}/disable_cdrom_repo.yml
- hosts: k8s-cluster:children
  become: true
  become_method: sudo
  tasks:
    - name: Disable the cdrom repository
      shell: sed -i '/deb cdrom:/s/^/#/g' /etc/apt/sources.list
    - name: Update cache
      shell: apt-get -y update; apt-get -y upgrade; apt-get -y autoremove
```
## Install Kubernetes
#### Clone Kubespray
```
cd ~
rm -rf kubespray
git clone https://github.com/kubernetes-incubator/kubespray.git
cd kubespray
git reset --hard e23fd5c
sudo pip install --upgrade pip
sudo pip install -r requirements.txt
cp -r inventory my_inventory
ln -sfn ${AWB_HOSTS} my_inventory/local/inventory
```
#### Edit `all.yml`
```
nano my_inventory/local/group_vars/all.yml
```
```
bootstrap_os: ubuntu
kubelet_load_modules: true
```
#### Edit `k8s-cluster.yml`
```
nano my_inventory/local/group_vars/k8s-cluster.yml
```
```
kube_version: v1.9.5
kube_network_plugin: calico
```
#### Edit `verify-settings.yml`
```
nano roles/kubernetes/preinstall/tasks/verify-settings.yml
```
```
- name: Stop if memory is too small for masters
  assert:
    that: ansible_memtotal_mb >= 500
  ignore_errors: "{{ ignore_assert_errors }}"
  when: inventory_hostname in groups['kube-master']
- name: Stop if memory is too small for nodes
  assert:
    that: ansible_memtotal_mb >= 500
  ignore_errors: "{{ ignore_assert_errors }}"
  when: inventory_hostname in groups['kube-node']
```
#### Edit `main.yml`
```
nano roles/docker/defaults/main.yml
```
```
docker_version: '17.03'
```
#### Edit main script
```
nano ~/kubespray/cluster.yml
```
```
- hosts: k8s-cluster:etcd:calico-rr
  any_errors_fatal: "{{ any_errors_fatal | default(true) }}"
  vars:
    - docker_dns_servers_strict: no
  roles:
    - { role: kubespray-defaults}
    - { role: kubernetes/preinstall, tags: preinstall }
    - { role: docker, tags: docker, when: manage_docker|default(true) }
    - role: rkt
      tags: rkt
      when: "'rkt' in [etcd_deployment_type, kubelet_deployment_type, vault_deployment_type]"
    - { role: download, tags: download, skip_downloads: false }
  environment: "{{proxy_env}}"
```
#### Install python packages
```
sudo pip install ansible netaddr
```
#### Manage `iptables`
```
sudo nano /etc/network/if-up.d/00-iptables
```
```
#!/bin/sh
iptables-restore < /etc/firewall.conf
ip6tables-restore < /etc/firewall6.conf
```
```
sudo chmod +x /etc/network/if-up.d/00-iptables
sudo iptables -A INPUT -p tcp --dport 2379 -j ACCEPT
sudo iptables -A INPUT -p tcp --dport 2380 -j ACCEPT
sudo iptables-save | sudo tee /etc/firewall.conf
sudo ip6tables-save | sudo tee /etc/firewall6.conf
```
#### Reboot & make cleanup
```
sudo rm -rf /tmp/*
sudo rm -rf /var/lib/etcd
sudo rm -rf /etc/systemd/system/etcd
sudo shutdown -r now

sudo swapoff -a
```
#### Deploy Kubespray
```
cd ~/kubespray
rm -f ~/kubespray.log
ansible-playbook --flush-cache -i ~/kubespray/my_inventory/local/inventory --ask-become-pass ~/kubespray/cluster.yml -b -vvv | tee -a ~/kubespray.log
```
#### Is Kubespray deployed?
```
kubectl get all
```
