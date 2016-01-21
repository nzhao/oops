let g:alternateSearchPath = 'wdr:include/spin,wdr:include/kron,wdr:include/quantum,dr:include/misc,wdr:source/spin/Spin,wdr:source/spin/SpinCollection,wdr:source/spin/SpinCluster,wdr:source/spin/SpinInteraction,wdr:source/quantum/QuantumOperator,wdr:source/quantum/QuantumState,wdr:source/quantum/QuantumEvolution,wdr:source/kron,wdr:source/misc'
command WM wa | make!
command RUN !../bin/cce
command DOC wa | make! doc
command TT tabe. | copen | wincmd k
