let g:alternateSearchPath = 'wdr:include/spin,wdr:include/kron,wdr:include/misc,wdr:source/spin/Spin,wdr:source/spin/SpinCollection,wdr:source/spin/SpinCluster,wdr:source/kron'
command WM wa | make 
command RUN !../bin/cce
command WMR wa | make | !../bin/cce 
