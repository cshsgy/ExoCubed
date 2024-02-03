#!/bin/bash

# The tarball and the destination directory
EIGEN="eigen-3.4.0"
TARFILE="$EIGEN.tar.gz"
DESTDIR="$HOME/opt/"

# Check if the destination directory exists, if not create it
if [ ! -d "$DESTDIR" ]; then
    mkdir -p "$DESTDIR"
fi

# Untar the tarball
tar -xvzf "$TARFILE"

# copy to destination
cp -r "$EIGEN" "$DESTDIR/eigen3"

# Print out a message when done
if [ $? -eq 0 ]; then
    echo "$EIGEN has been successfully installed to $DESTDIR"
else
    echo "There was a problem during the installation."
fi
