/*
 * circuit.h - the data structures that describe a circuit.
 *
 * A netlist line like "R1 in out 1k" becomes one Element. Node names are
 * turned into integers: node 0 is ground, nodes 1..n-1 are everything else.
 */
#ifndef CIRCUIT_H
#define CIRCUIT_H

#define MAX_NAME 32

typedef enum {
    EL_R,   /* resistor                          R n+ n- value            */
    EL_C,   /* capacitor                         C n+ n- value            */
    EL_L,   /* inductor                          L n+ n- value            */
    EL_V,   /* independent voltage source        V n+ n- <source spec>    */
    EL_I,   /* independent current source        I n+ n- <source spec>    */
    EL_D,   /* diode (anode, cathode)            D a  k  model            */
    EL_S,   /* voltage-controlled switch         S n+ n- nc+ nc- model    */
    EL_E    /* voltage-controlled voltage source E n+ n- nc+ nc- gain     */
} ElemType;

typedef enum { SRC_DC, SRC_PULSE, SRC_SIN, SRC_PWL } SrcKind;

/* Waveform of an independent source. */
typedef struct {
    SrcKind kind;
    double dc;          /* DC value                                         */
    double p[7];        /* PULSE: v1 v2 td tr tf pw per | SIN: vo va f td    */
    int npwl;           /* PWL: number of (time, value) pairs                */
    double *pwl_t;
    double *pwl_v;
} Source;

typedef enum { MOD_D, MOD_SW } ModelType;

typedef struct {
    char name[MAX_NAME];
    ModelType type;
    double is, n;           /* diode: saturation current, ideality factor   */
    double ron, roff, vt;   /* switch: on/off resistance, threshold voltage */
} Model;

typedef struct {
    ElemType type;
    char name[MAX_NAME];
    int n1, n2;         /* the two terminals                                */
    int nc1, nc2;       /* controlling nodes (S and E only)                 */
    double value;       /* R: ohms, C: farads, L: henries, E: gain          */
    Source src;         /* V and I only                                     */
    int model;          /* index into Circuit.models (D and S)              */
    int branch;         /* extra unknown (branch current) for V, L, E; -1 if none */

    /* integration state, updated after every accepted time step */
    double v_prev;      /* voltage across the element at the last step     */
    double i_prev;      /* current through the element at the last step    */
    double vd_iter;     /* diode: voltage used in the previous Newton iteration */
} Element;

typedef struct {
    int nnodes;                 /* including ground                        */
    char **node_names;
    Element *el;
    int nel, cap_el;
    Model *models;
    int nmodels, cap_models;
    int nbranch;                /* number of branch-current unknowns       */
    /* .tran tstep tstop [tstart] */
    int has_tran;
    double tstep, tstop, tstart;
    char title[256];
} Circuit;

/* number of MNA unknowns: node voltages (minus ground) + branch currents */
static inline int circuit_size(const Circuit *c) { return c->nnodes - 1 + c->nbranch; }

void circuit_free(Circuit *c);

#endif
