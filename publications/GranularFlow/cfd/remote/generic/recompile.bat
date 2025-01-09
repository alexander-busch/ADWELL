if exist "libudf\win64\2ddp\makefile" cd libudf\win64\2ddp & nmake clean & nmake & cd ..\..\.. &
if exist "libudf\win64\2ddp_host\makefile" cd libudf\win64\2ddp_host & nmake clean & nmake & cd ..\..\.. & cd libudf\win64\2ddp_node & nmake clean & nmake & cd ..\..\.. &
if exist "libudf\win64\3ddp\makefile" cd libudf\win64\3ddp & nmake clean & nmake & cd ..\..\.. &
if exist "libudf\win64\3ddp_host\makefile" cd libudf\win64\3ddp_host & nmake clean & nmake & cd ..\..\.. & cd libudf\win64\3ddp_node & nmake clean & nmake & cd ..\..\.. &