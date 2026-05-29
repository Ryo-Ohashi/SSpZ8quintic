This repository provides the **Magma** implementation of the algorithm proposed in the following paper:

> [Superspecial plane quintics with large automorphism groups](https://arxiv.org/abs/2605.29624) by Ryo Ohashi

This code was used to generate the data presented in Table 3 of the paper. Running

```
load "main.m";
```

outputs, for each prime $p$ satisfying $13 < p < 10000$ and $p \equiv 7 \pmod{8}$, the number of isomorphism classes of superspecial plane quintic curves with a cyclic automorphism group of order 8 in characteristic $p$.
