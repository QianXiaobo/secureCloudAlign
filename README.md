Secure Cloud Align

The elastic and inexpensive computing resources such as clouds have been recognized as a useful solution to analyzing massive human genomic data (e.g., acquired by using next-generation sequencers) in biomedical research. However, outsourcing human genome computation to public or commercial clouds were hindered due to privacy concerns: even a small number of human genome sequences contain sufficient information for identifying the donor of the genomic data. This issue cannot be directly addressed by existing security and cryptographic techniques (such as homomorphic encryption), because they are too heavyweight to carry out practical genome computation tasks on massive data. In this paper, we present a secure algorithm to accomplish the read mapping, one of the most basic tasks in human genomic data analysis based on a {\em hybrid cloud computing} model. Comparing with the existing approaches, our algorithm delegates most computation to the public cloud, while only performing encryption and decryption on the private cloud, and thus makes the maximum use of the computing resource of the public cloud. Furthermore, our algorithm reports similar results as the non-secure reads mapping algorithms, including the alignment between reads and the reference genome, which can be directly used in the downstream analysis such as the inference of genomic variations. We implemented the algorithm in C++ and Python on a hybrid cloud system, in which the public cloud uses an Apache Spark system. The implementation is released as open source software at github (https://github.com/zhaoyanswill/secureCloudAlign).
