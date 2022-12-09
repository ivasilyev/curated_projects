# Start server

```shell script
export IMG="redis:latest" && \
docker pull "${IMG}" && \
docker run \
    --expose 6379/tcp \
    --interactive \
    --net=host \
    --rm \
    --tty \
    "${IMG}" \
    redis-server
```

# Check if the server works (from separate console)

```shell script
echo PING | nc localhost 6379 -w 1
```

```text
+PONG
```

# Start client

```shell script
export IMG="redis:latest" && \
docker pull "${IMG}" && \
docker run \
    --interactive \
    --net=host \
    --rm \
    --tty \
    "${IMG}" \
    bash
```

# Refill database

```shell script
alias REDIS="redis-cli -h localhost -p 6379"

REDIS flushall
echo
REDIS rpush test_list '{"key1": "value1"}'
REDIS rpush test_list '{"key2": "value2"}'
echo
REDIS lrange test_list 0 -1
```

```text
OK

(integer) 1
(integer) 2

1) "{\"key1\": \"value1\"}"
2) "{\"key2\": \"value2\"}"
```

# Blocking move from left to right

```shell script
# BLMOVE source destination <LEFT | RIGHT> <LEFT | RIGHT> timeout
REDIS blmove test_list test_list:processing left right 2
echo
REDIS lrange test_list 0 -1
echo
REDIS lrange test_list:processing 0 -1
```

```text
"{\"key1\": \"value1\"}"

1) "{\"key2\": \"value2\"}"

1) "{\"key1\": \"value1\"}"
```

# Blocking move from right to left

```shell script
REDIS blmove test_list:processing test_list right left 2
echo
REDIS lrange test_list 0 -1
echo
REDIS lrange test_list:processing 0 -1
```

```text
"{\"key1\": \"value1\"}"

1) "{\"key1\": \"value1\"}"
2) "{\"key2\": \"value2\"}"

(empty array)
```
