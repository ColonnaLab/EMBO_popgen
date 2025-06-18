# SERVER CONNECTION 2025

The Server can be accessed by speakers and students at:
```
embo-popgen25.cloud.ba.infn.it
```

The software has been installed system-wide and should be visible by everyone. 
To access the server, please follow one of the two below options:

1. via SSH:
using your own username and private key (the one corresponding to the public key which you provided us with); SSH with password is not allowed; using this command

````
ssh -i <PATH_TO_PRIVATE_KEY> embo-popgen25.cloud.ba.infn.it -l <USERNAME>
````

(the "-i <PATH_TO_PRIVATE_KEY>" flag can be omitted if the private key is in the standard location, in the $HOME/.ssh/ directory on UNIX-like systems as Linux, Mac or WSL in Windows)

2. via WEB:
 using the RStudio Browser interface only (which is just an  authenticated RStudio web interface), pointing the   browser to the URL:

   https://embo-popgen25.cloud.ba.infn.it/

The usernames are the same used for SSH, but in this case you must insert your own password to access the system.

**Usernames and passwords will be provided via Slack: if you have not received the invitation to join the EMBO POPGEN 2025 channel please get in touch with us ASAP.**
