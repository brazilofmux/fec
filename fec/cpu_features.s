.section .text
.global cpu_features
        .type cpu_features,@function
cpu_features:
        pushq %rbx
        pushq %rcx
        pushq %rdx
        movl $1, %eax
        cpuid
        movl %edx, %eax
        popq %rdx
        popq %rcx
        popq %rbx
        ret
        .size cpu_features, .-cpu_features
        .section .note.GNU-stack,"",@progbits
