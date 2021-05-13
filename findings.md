# T=1

## Unrooted counterexample

### 007

Turns a speciation into a lower duplication, inverting an existing transfer and saving two segmental losses.

THL-SPFS = 3 (1 + 2)

`SuperReconciliation(input=SuperReconciliationInput(synteny_tree="(0_0,(1_1,0_2)1_10_2)0_01_10_2;" species_tree="(0,1)01;" leaf_labeling="0_0:a,1_1:a,0_2:ab") reconciliation="0_01_10_2:01,0_0:0,1_10_2:1,1_1:1,0_2:0" labeling="0_01_10_2:ab,0_0:a,1_10_2:ab,1_1:a,0_2:ab")`

BF = 2 (2 + 0)

`SuperReconciliation(input=SuperReconciliationInput(synteny_tree="(0_0,(1_1,0_2)1_10_2)0_01_10_2;" species_tree="(0,1)01;" leaf_labeling="0_0:a,1_1:a,0_2:ab") reconciliation="0_01_10_2:0,0_0:0,1_10_2:0,1_1:1,0_2:0" labeling="0_01_10_2:ab,0_0:a,1_10_2:ab,1_1:a,0_2:ab")`

### 010

Turns a speciation into a lower transfer, saving two segmental losses.

THL-SPFS = 5 (1 + 4)

`SuperReconciliation(input=SuperReconciliationInput(synteny_tree="((0_0,1_1)0_01_1,(0_2,1_3)0_21_3)0_01_10_21_3;" species_tree="(0,1)01;" leaf_labeling="0_0:a,1_1:b,0_2:a,1_3:b") reconciliation="0_01_10_21_3:01,0_01_1:01,0_0:0,1_1:1,0_21_3:01,0_2:0,1_3:1" labeling="0_01_10_21_3:ab,0_01_1:ab,0_0:a,1_1:b,0_21_3:ab,0_2:a,1_3:b")`

BF = 4 (2 + 2)

`SuperReconciliation(input=SuperReconciliationInput(synteny_tree="((0_0,1_1)0_01_1,(0_2,1_3)0_21_3)0_01_10_21_3;" species_tree="(0,1)01;" leaf_labeling="0_0:a,1_1:b,0_2:a,1_3:b") reconciliation="0_01_10_21_3:01,0_01_1:0,0_0:0,1_1:1,0_21_3:1,0_2:0,1_3:1" labeling="0_01_10_21_3:ab,0_01_1:ab,0_0:a,1_1:b,0_21_3:ab,0_2:a,1_3:b")`

## Rooted counterexample

### 008

Turns two speciations into lower transfers, saving segmental losses.

THL-SPFS = 5 (1 + 4)

`SuperReconciliation(input=SuperReconciliationInput(synteny_tree="((0_0,1_1)0_01_1,(0_2,1_3)0_21_3)0_01_10_21_3;" species_tree="(0,1)01;" leaf_labeling="0_0:a,1_1:b,0_2:a,1_3:b") reconciliation="0_01_10_21_3:01,0_01_1:01,0_0:0,1_1:1,0_21_3:01,0_2:0,1_3:1" labeling="0_01_10_21_3:ab,0_01_1:ab,0_0:a,1_1:b,0_21_3:ab,0_2:a,1_3:b")`

BF = 4 (2 + 2)

`SuperReconciliation(input=SuperReconciliationInput(synteny_tree="((0_0,1_1)0_01_1,(0_2,1_3)0_21_3)0_01_10_21_3;" species_tree="(0,1)01;" leaf_labeling="0_0:a,1_1:b,0_2:a,1_3:b") reconciliation="0_01_10_21_3:01,0_01_1:0,0_0:0,1_1:1,0_21_3:1,0_2:0,1_3:1" labeling="0_01_10_21_3:ab,0_01_1:ab,0_0:a,1_1:b,0_21_3:ab,0_2:a,1_3:b")`

### 009

Turns a root speciation into a higher duplication at the cost of an added full loss. Allows for turning a transfer back into a speciation and for less segmental losses.

THL-SPFS = 3 (1 + 2)

`SuperReconciliation(input=SuperReconciliationInput(synteny_tree="((0_0,1_1)0_01_1,0_2)0_01_10_2;" species_tree="(0,1)01;" leaf_labeling="0_0:b,1_1:b,0_2:abc") reconciliation="0_01_10_2:01,0_01_1:1,0_0:0,1_1:1,0_2:0" labeling="0_01_10_2:abc,0_01_1:b,0_0:b,1_1:b,0_2:abc")`

BF = 2 (2 + 0)

`SuperReconciliation(input=SuperReconciliationInput(synteny_tree="((0_0,1_1)0_01_1,0_2)0_01_10_2;" species_tree="(0,1)01;" leaf_labeling="0_0:b,1_1:b,0_2:abc") reconciliation="0_01_10_2:01,0_01_1:01,0_0:0,1_1:1,0_2:0" labeling="0_01_10_2:abc,0_01_1:b,0_0:b,1_1:b,0_2:abc")`

# T=2

### 011

THL-SPFS = 4 (2 + 2)

`SuperReconciliation(input=SuperReconciliationInput(synteny_tree="(0_0,(1_1,0_2)1_10_2)0_01_10_2;", species_tree="(0,1)01;", leaf_labeling="0_0:a,1_1:b,0_2:abc"), reconciliation="0_01_10_2:01,0_0:0,1_10_2:01,1_1:1,0_2:0", labeling="0_01_10_2:abc,0_0:a,1_10_2:abc,1_1:b,0_2:abc")`

BF = 3 (3 + 0)

`SuperReconciliation(input=SuperReconciliationInput(synteny_tree="(0_0,(1_1,0_2)1_10_2)0_01_10_2;", species_tree="(0,1)01;", leaf_labeling="0_0:a,1_1:b,0_2:abc"), reconciliation="0_01_10_2:0,0_0:0,1_10_2:0,1_1:1,0_2:0", labeling="0_01_10_2:abc,0_0:a,1_10_2:abc,1_1:b,0_2:abc")`

# T=3

### 012

THL-SPFS = 5 (2 + 3)

`SuperReconciliation(input=SuperReconciliationInput(synteny_tree="((0_0,1_1)0_01_1,0_2)0_01_10_2;", species_tree="(0,1)01;", leaf_labeling="0_0:b,1_1:ac,0_2:abc"), reconciliation="0_01_10_2:01,0_01_1:01,0_0:0,1_1:1,0_2:0", labeling="0_01_10_2:abc,0_01_1:abc,0_0:b,1_1:ac,0_2:abc")`

BF = 4 (3 + 1)

`SuperReconciliation(input=SuperReconciliationInput(synteny_tree="((0_0,1_1)0_01_1,0_2)0_01_10_2;", species_tree="(0,1)01;", leaf_labeling="0_0:b,1_1:ac,0_2:abc"), reconciliation="0_01_10_2:01,0_01_1:1,0_0:0,1_1:1,0_2:0", labeling="0_01_10_2:abc,0_01_1:abc,0_0:b,1_1:ac,0_2:abc")`

# T=4

??

### 013

THL-SPFS = 7 (4 + 3)

`SuperReconciliation(input=SuperReconciliationInput(synteny_tree="((0_0,1_1)0_01_1,1_2)0_01_11_2;", species_tree="(0,(1,2)12)012;", leaf_labeling="0_0:a,1_1:bc,1_2:ac"), reconciliation="0_01_11_2:012,0_01_1:012,0_0:0,1_1:1,1_2:1", labeling="0_01_11_2:abc,0_01_1:abc,0_0:a,1_1:bc,1_2:ac")`

BF = 6 (5 + 1)

`SuperReconciliation(input=SuperReconciliationInput(synteny_tree="((0_0,1_1)0_01_1,1_2)0_01_11_2;", species_tree="(0,(1,2)12)012;", leaf_labeling="0_0:a,1_1:bc,1_2:ac"), reconciliation="0_01_11_2:1,0_01_1:1,0_0:0,1_1:1,1_2:1", labeling="0_01_11_2:bac,0_01_1:bac,0_0:a,1_1:bc,1_2:ac")`

