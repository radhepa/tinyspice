/*
 * netlist.c - SPICE-style netlist parser.
 *
 * Supported:
 *   first line            title
 *   * ...                 comment line;  ";" starts an end-of-line comment
 *   + ...                 continues the previous line
 *   Rxx/Cxx/Lxx n+ n- value
 *   Vxx/Ixx n+ n- [dc] value | pulse(v1 v2 td tr tf pw per) | sin(vo va f [td]) | pwl(t1 v1 t2 v2 ...)
 *   Dxx anode cathode model
 *   Sxx n+ n- nc+ nc- model
 *   Exx n+ n- nc+ nc- gain
 *   .model name d(is=... n=...)  |  .model name sw(ron=... roff=... vt=...)
 *   .tran tstep tstop [tstart]
 *   .end
 */
#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "netlist.h"

#define MAX_LINE 4096
#define MAX_TOK 256

/* ----------------------------------------------------------------- helpers */

static void set_err(char *err, size_t n, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(err, n, fmt, ap);
    va_end(ap);
}

int parse_value(const char *s, double *out)
{
    char *end;
    double v = strtod(s, &end);
    if (end == s)
        return -1;
    /* SPICE suffixes are case-insensitive; "meg" must be checked before "m" */
    char suf[8] = {0};
    for (int i = 0; i < 7 && end[i]; i++)
        suf[i] = (char)tolower((unsigned char)end[i]);
    if (strncmp(suf, "meg", 3) == 0)      v *= 1e6;
    else if (suf[0] == 'f')               v *= 1e-15;
    else if (suf[0] == 'p')               v *= 1e-12;
    else if (suf[0] == 'n')               v *= 1e-9;
    else if (suf[0] == 'u')               v *= 1e-6;
    else if (suf[0] == 'm')               v *= 1e-3;
    else if (suf[0] == 'k')               v *= 1e3;
    else if (suf[0] == 'g')               v *= 1e9;
    else if (suf[0] == 't')               v *= 1e12;
    /* anything after the suffix (units such as "F" or "ohm") is ignored */
    *out = v;
    return 0;
}

/* Split a line into tokens. Parentheses and commas are separators, and '='
 * becomes its own token, so "d(is=1e-14)" -> "d" "is" "=" "1e-14".
 * `buf` receives a spaced-out copy of the line; tokens point into it. */
static int tokenize(const char *line, char *buf, size_t buflen, char *tok[], int max)
{
    size_t k = 0;
    for (const char *p = line; *p && k + 4 < buflen; p++) {
        if (*p == '(' || *p == ')' || *p == ',') {
            buf[k++] = ' ';
        } else if (*p == '=') {
            buf[k++] = ' ';
            buf[k++] = '=';
            buf[k++] = ' ';
        } else {
            buf[k++] = *p;
        }
    }
    buf[k] = '\0';

    int n = 0;
    for (char *t = strtok(buf, " \t"); t && n < max; t = strtok(NULL, " \t"))
        tok[n++] = t;
    return n;
}

static char *str_dup(const char *s)
{
    size_t n = strlen(s) + 1;
    char *p = malloc(n);
    memcpy(p, s, n);
    return p;
}

static void lower(char *s)
{
    for (; *s; s++)
        *s = (char)tolower((unsigned char)*s);
}

static int node_index(Circuit *c, const char *name)
{
    if (strcmp(name, "0") == 0 || strcmp(name, "gnd") == 0)
        return 0;
    for (int i = 1; i < c->nnodes; i++)
        if (strcmp(c->node_names[i], name) == 0)
            return i;
    c->node_names = realloc(c->node_names, sizeof(char *) * (size_t)(c->nnodes + 1));
    c->node_names[c->nnodes] = str_dup(name);
    return c->nnodes++;
}

static Element *new_element(Circuit *c)
{
    if (c->nel == c->cap_el) {
        c->cap_el = c->cap_el ? 2 * c->cap_el : 32;
        c->el = realloc(c->el, sizeof(Element) * (size_t)c->cap_el);
    }
    Element *e = &c->el[c->nel++];
    memset(e, 0, sizeof *e);
    e->branch = -1;
    e->model = -1;
    return e;
}

static Model *new_model(Circuit *c)
{
    if (c->nmodels == c->cap_models) {
        c->cap_models = c->cap_models ? 2 * c->cap_models : 8;
        c->models = realloc(c->models, sizeof(Model) * (size_t)c->cap_models);
    }
    Model *m = &c->models[c->nmodels++];
    memset(m, 0, sizeof *m);
    return m;
}

/* ------------------------------------------------------- source waveforms */

static int parse_source(char **tok, int n, Source *s, char *err, size_t errlen)
{
    memset(s, 0, sizeof *s);
    s->kind = SRC_DC;
    int i = 0;
    while (i < n) {
        if (strcmp(tok[i], "dc") == 0 && i + 1 < n) {
            if (parse_value(tok[i + 1], &s->dc)) goto bad;
            i += 2;
        } else if (strcmp(tok[i], "pulse") == 0) {
            s->kind = SRC_PULSE;
            int k;
            for (k = 0; k < 7 && i + 1 + k < n; k++)
                if (parse_value(tok[i + 1 + k], &s->p[k])) break;
            if (k < 2) goto bad;
            s->dc = s->p[0];
            i += 1 + k;
        } else if (strcmp(tok[i], "sin") == 0) {
            s->kind = SRC_SIN;
            int k;
            for (k = 0; k < 4 && i + 1 + k < n; k++)
                if (parse_value(tok[i + 1 + k], &s->p[k])) break;
            if (k < 3) goto bad;
            s->dc = s->p[0];
            i += 1 + k;
        } else if (strcmp(tok[i], "pwl") == 0) {
            s->kind = SRC_PWL;
            int k = i + 1;
            s->pwl_t = malloc(sizeof(double) * (size_t)n);
            s->pwl_v = malloc(sizeof(double) * (size_t)n);
            while (k + 1 < n && !parse_value(tok[k], &s->pwl_t[s->npwl])
                   && !parse_value(tok[k + 1], &s->pwl_v[s->npwl])) {
                s->npwl++;
                k += 2;
            }
            if (s->npwl < 1) goto bad;
            s->dc = s->pwl_v[0];
            i = k;
        } else if (!parse_value(tok[i], &s->dc)) {
            i++;                              /* bare value = DC */
        } else {
            goto bad;
        }
    }
    return 0;
bad:
    set_err(err, errlen, "cannot read source specification near '%s'", i < n ? tok[i] : "(end)");
    return -1;
}

/* ---------------------------------------------------------------- parser */

static int parse_model(Circuit *c, char **tok, int n, char *err, size_t errlen)
{
    if (n < 3) {
        set_err(err, errlen, ".model needs a name and a type");
        return -1;
    }
    Model *m = new_model(c);
    snprintf(m->name, MAX_NAME, "%s", tok[1]);
    if (strcmp(tok[2], "d") == 0) {
        m->type = MOD_D;
        m->is = 1e-14;
        m->n = 1.0;
    } else if (strcmp(tok[2], "sw") == 0) {
        m->type = MOD_SW;
        m->ron = 1.0;
        m->roff = 1e12;
        m->vt = 0.0;
    } else {
        set_err(err, errlen, "unknown model type '%s' (use d or sw)", tok[2]);
        return -1;
    }
    for (int i = 3; i + 2 < n; i += 3) {       /* key = value triples */
        if (strcmp(tok[i + 1], "=") != 0)
            break;
        double v;
        if (parse_value(tok[i + 2], &v)) {
            set_err(err, errlen, "bad value in .model %s: %s", m->name, tok[i + 2]);
            return -1;
        }
        const char *k = tok[i];
        if (strcmp(k, "is") == 0) m->is = v;
        else if (strcmp(k, "n") == 0) m->n = v;
        else if (strcmp(k, "ron") == 0) m->ron = v;
        else if (strcmp(k, "roff") == 0) m->roff = v;
        else if (strcmp(k, "vt") == 0) m->vt = v;
        /* other SPICE parameters are accepted and ignored */
    }
    return 0;
}

static int parse_line(Circuit *c, char *line, char **model_refs, int lineno, char *err, size_t errlen)
{
    char *tok[MAX_TOK];
    static char buf[MAX_LINE * 8];
    lower(line);
    int n = tokenize(line, buf, sizeof buf, tok, MAX_TOK);
    if (n == 0)
        return 0;

    if (tok[0][0] == '.') {
        if (strcmp(tok[0], ".model") == 0)
            return parse_model(c, tok, n, err, errlen);
        if (strcmp(tok[0], ".tran") == 0) {
            if (n < 3 || parse_value(tok[1], &c->tstep) || parse_value(tok[2], &c->tstop)) {
                set_err(err, errlen, "line %d: .tran tstep tstop [tstart]", lineno);
                return -1;
            }
            if (n > 3) parse_value(tok[3], &c->tstart);
            c->has_tran = 1;
        }
        return 0;   /* .op, .end, .options ... are accepted silently */
    }

    Element *e = new_element(c);
    snprintf(e->name, MAX_NAME, "%s", tok[0]);
    int need = (tok[0][0] == 's' || tok[0][0] == 'e') ? 6 : 4;
    if (n < need) {
        set_err(err, errlen, "line %d: '%s' needs at least %d fields", lineno, tok[0], need);
        return -1;
    }
    e->n1 = node_index(c, tok[1]);
    e->n2 = node_index(c, tok[2]);

    switch (tok[0][0]) {
    case 'r': case 'c': case 'l':
        e->type = tok[0][0] == 'r' ? EL_R : tok[0][0] == 'c' ? EL_C : EL_L;
        if (parse_value(tok[3], &e->value) || (e->type == EL_R && e->value == 0.0)) {
            set_err(err, errlen, "line %d: bad value '%s'", lineno, tok[3]);
            return -1;
        }
        if (e->type == EL_L)
            e->branch = c->nbranch++;
        break;
    case 'v': case 'i':
        e->type = tok[0][0] == 'v' ? EL_V : EL_I;
        if (parse_source(tok + 3, n - 3, &e->src, err, errlen))
            return -1;
        if (e->type == EL_V)
            e->branch = c->nbranch++;
        break;
    case 'd':
        e->type = EL_D;
        model_refs[c->nel - 1] = str_dup(tok[3]);
        break;
    case 's':
        e->type = EL_S;
        e->nc1 = node_index(c, tok[3]);
        e->nc2 = node_index(c, tok[4]);
        model_refs[c->nel - 1] = str_dup(tok[5]);
        break;
    case 'e':
        e->type = EL_E;
        e->nc1 = node_index(c, tok[3]);
        e->nc2 = node_index(c, tok[4]);
        if (parse_value(tok[5], &e->value)) {
            set_err(err, errlen, "line %d: bad gain '%s'", lineno, tok[5]);
            return -1;
        }
        e->branch = c->nbranch++;
        break;
    default:
        set_err(err, errlen, "line %d: unsupported element '%s'", lineno, tok[0]);
        return -1;
    }
    return 0;
}

int parse_netlist(const char *path, Circuit *c, char *err, size_t errlen)
{
    FILE *f = fopen(path, "r");
    if (!f) {
        set_err(err, errlen, "cannot open '%s'", path);
        return -1;
    }
    memset(c, 0, sizeof *c);
    c->nnodes = 1;
    c->node_names = malloc(sizeof(char *));
    c->node_names[0] = str_dup("0");

    /* read the whole file, joining "+" continuation lines */
    char raw[MAX_LINE], logical[MAX_LINE * 4] = "";
    int lineno = 0, first = 1, logical_line = 0, rc = 0;
    int cap_refs = 1024;
    char **model_refs = calloc((size_t)cap_refs, sizeof(char *));

    for (;;) {
        char *got = fgets(raw, sizeof raw, f);
        if (got) {
            lineno++;
            raw[strcspn(raw, "\r\n")] = '\0';
            char *semi = strchr(raw, ';');
            if (semi) *semi = '\0';
            if (first) {                      /* SPICE: line 1 is the title */
                snprintf(c->title, sizeof c->title, "%.255s", raw);
                first = 0;
                continue;
            }
            if (raw[0] == '+') {              /* continuation */
                strncat(logical, " ", sizeof logical - strlen(logical) - 1);
                strncat(logical, raw + 1, sizeof logical - strlen(logical) - 1);
                continue;
            }
        }
        /* a new line starts (or EOF): process the previous logical line */
        if (logical[0] && logical[0] != '*') {
            if (c->nel + 1 >= cap_refs) {
                model_refs = realloc(model_refs, sizeof(char *) * (size_t)(cap_refs * 2));
                memset(model_refs + cap_refs, 0, sizeof(char *) * (size_t)cap_refs);
                cap_refs *= 2;
            }
            char tmp[MAX_LINE * 4];
            snprintf(tmp, sizeof tmp, "%s", logical);
            char *p = tmp;
            while (isspace((unsigned char)*p)) p++;
            if (strncmp(p, ".end", 4) == 0 && (p[4] == '\0' || isspace((unsigned char)p[4]))) {
                logical[0] = '\0';
                break;
            }
            if (parse_line(c, tmp, model_refs, logical_line, err, errlen)) {
                rc = -1;
                break;
            }
        }
        if (!got)
            break;
        snprintf(logical, sizeof logical, "%s", raw);
        logical_line = lineno;
    }
    fclose(f);

    /* resolve model names now that every .model line has been read */
    for (int i = 0; rc == 0 && i < c->nel; i++) {
        if (!model_refs[i])
            continue;
        Element *e = &c->el[i];
        for (int m = 0; m < c->nmodels; m++)
            if (strcmp(c->models[m].name, model_refs[i]) == 0)
                e->model = m;
        if (e->model < 0) {
            set_err(err, errlen, "%s: model '%s' not defined", e->name, model_refs[i]);
            rc = -1;
        } else if ((e->type == EL_D) != (c->models[e->model].type == MOD_D)) {
            set_err(err, errlen, "%s: model '%s' has the wrong type", e->name, model_refs[i]);
            rc = -1;
        }
    }
    for (int i = 0; i < cap_refs; i++)
        free(model_refs[i]);
    free(model_refs);
    return rc;
}

void circuit_free(Circuit *c)
{
    for (int i = 0; i < c->nnodes; i++)
        free(c->node_names[i]);
    free(c->node_names);
    for (int i = 0; i < c->nel; i++) {
        free(c->el[i].src.pwl_t);
        free(c->el[i].src.pwl_v);
    }
    free(c->el);
    free(c->models);
    memset(c, 0, sizeof *c);
}
