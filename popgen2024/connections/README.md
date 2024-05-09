# SERVER CONNECTION

The Server can be accessed by speakers and students at: 

```
embo-popgen-2024.recas.ba.infn.it 
```
The software has been installed system-wide and should be visible by everyone.

To access the server, please follow one of the two below options:

1. via SSH: 
using your own username and private key (the one corresponding to the public key which you provided us with); SSH with password is not allowed; the command to be used is

```
ssh -i <PATH_TO_PRIVATE_KEY> embo-popgen-2024.recas.ba.infn.it -l <USERNAME>
```

(the "-i <PATH_TO_PRIVATE_KEY>" flag can be omitted if the private key
is the standard location $HOME/.ssh/id_rsa on UNIX-like systems as
Linux, Mac or WSL in Windows)

2. via WEB:
using the JupyterHub interface only (which is just an authenticated jupyterlab notebook web interface), pointing the browser to the URL: https://embo-popgen-2023.recas.ba.infn.it

The usernames are the same used for SSH, but in this case you must insert your own password to access the system.

**Usernames and passwords will be provided via Slack: if you have not received the invitation to join the EMBO POPGEN 2023 channel please get in touch with us ASAP.**
