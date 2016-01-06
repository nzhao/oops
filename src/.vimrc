let g:alternateSearchPath = 'wdr:include/spin,wdr:include/kron,wdr:include/misc,wdr:source/spin/Spin,wdr:source/spin/SpinCollection,wdr:source/spin/SpinCluster,wdr:source/spin/SpinInteraction,wdr:source/kron'
command WM wa | make 
command RUN !../bin/cce
command DOC wa | make doc
