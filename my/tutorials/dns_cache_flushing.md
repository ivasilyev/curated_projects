# DNS Cache Flushing

* Connect to the network that does not have a poisoned DNS cache;
* Use the commands:

```
# Windows

ipconfig /flushdns

# Linux
sudo systemd-resolve --flush-caches
sudo resolvectl flush-caches
```
