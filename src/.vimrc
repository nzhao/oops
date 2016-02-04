let g:alternateSearchPath = '
            \wdr:include/kron,
            \wdr:include/misc,
            \wdr:include/quantum,
            \wdr:include/spin,
            \wdr:source/kron,
            \wdr:source/misc,
            \wdr:source/quantum/QuantumOperator,
            \wdr:source/quantum/QuantumState,
            \wdr:source/quantum/QuantumEvolution,
            \wdr:source/spin/Spin,
            \wdr:source/spin/SpinCollection,
            \wdr:source/spin/SpinCluster,
            \wdr:source/spin/SpinInteraction'
command WM wa | make!
command RUN !../bin/cce
command DOC wa | make! doc
command TT tabe. | copen | wincmd k
