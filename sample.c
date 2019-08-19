/*  sample.c -- group data by sample.

    Copyright (C) 2010, 2011 Broad Institute.
    Copyright (C) 2013 Genome Research Ltd.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <config.h>

#include <stdlib.h>
#include <string.h>
#include "sample.h"
#include "htslib/khash.h"
KHASH_MAP_INIT_STR(sm, int)

bam_sample_t *bam_smpl_init(void)
{
    bam_sample_t *s;
    s = calloc(1, sizeof(bam_sample_t));
    s->rg2smid = kh_init(sm);
    s->sm2id = kh_init(sm);
    return s;
}

void bam_smpl_destroy(bam_sample_t *sm)
{
    int i;
    khint_t k;
    khash_t(sm) *rg2smid = (khash_t(sm)*)sm->rg2smid;
    if (sm == 0) return;
    for (i = 0; i < sm->n; ++i) free(sm->smpl[i]);
    free(sm->smpl);
    for (k = kh_begin(rg2smid); k != kh_end(rg2smid); ++k)
        if (kh_exist(rg2smid, k)) free((char*)kh_key(rg2smid, k));
    kh_destroy(sm, sm->rg2smid);
    kh_destroy(sm, sm->sm2id);
    free(sm);
}

static int add_pair(bam_sample_t *sm, khash_t(sm) *sm2id, const char *key, const char *val)
{
    khint_t k_rg, k_sm;
    char *k = NULL;
    int ret;
    khash_t(sm) *rg2smid = (khash_t(sm)*)sm->rg2smid;
    k_rg = kh_get(sm, rg2smid, key);
    if (k_rg != kh_end(rg2smid)) return 0; // duplicated @RG-ID
    k = strdup(key);
    if (!k) return -1;
    k_rg = kh_put(sm, rg2smid, k, &ret);
    if (ret < 0) {
        free(k);
        return -1;
    }
    k_sm = kh_get(sm, sm2id, val);
    if (k_sm == kh_end(sm2id)) { // absent
        if (sm->n == sm->m) {
            int new_m = sm->m ? sm->m * 2 : 1;
            char **new_smpl = realloc(sm->smpl, sizeof(char*) * new_m);
            if (!new_smpl) return -1;
            sm->m = new_m;
            sm->smpl = new_smpl;
        }
        sm->smpl[sm->n] = strdup(val);
        if (!sm->smpl[sm->n]) return -1;
        k_sm = kh_put(sm, sm2id, sm->smpl[sm->n], &ret);
        if (ret < 0) {
            free(sm->smpl[sm->n]);
            return -1;
        }
        kh_val(sm2id, k_sm) = sm->n++;
    }
    kh_val(rg2smid, k_rg) = kh_val(sm2id, k_sm);
    return 0;
}

int bam_smpl_add(bam_sample_t *sm, const char *fn, sam_hdr_t *hdr)
{
    kstring_t buf = KS_INITIALIZE, first_sm = KS_INITIALIZE;
    kstring_t id_val = KS_INITIALIZE, sm_val = KS_INITIALIZE;
    int n = 0, nrg, i;
    khash_t(sm) *sm2id = (khash_t(sm)*)sm->sm2id;
    if (hdr == NULL) {
        return add_pair(sm, sm2id, fn, fn);
    }
    nrg = sam_hdr_count_lines(hdr, "RG");
    if (nrg < 0) return -1;

    for (i = 0; i < nrg; i++) {
        // This follows earlier code that would break on @RG with no ID/SM.
        // The header API should enforce ID, but SM can still be absent.
        // Possibly it should add a default value for such lines instead?
        if (sam_hdr_find_tag_pos(hdr, "RG", i, "ID", &id_val) < 0) break;
        if (sam_hdr_find_tag_pos(hdr, "RG", i, "SM", &sm_val) < 0) break;
        if (kputs(fn, ks_clear(&buf)) < 0) goto fail;
        if (kputc('/', &buf) < 0) goto fail;
        if (kputs(id_val.s, &buf) < 0) goto fail;
        if (add_pair(sm, sm2id, buf.s, sm_val.s) < 0) goto fail;
        if (!first_sm.s) {
            first_sm.m = sm_val.m;
            first_sm.l = sm_val.l;
            first_sm.s = ks_release(&sm_val);
        }
        ++n;
    }

    if (n == 0) {
        if (add_pair(sm, sm2id, fn, fn) < 0) goto fail;
    } else if ( n==1 && first_sm.s ) {
        // If there is only one RG tag present in the header and reads are not
        // annotated, don't refuse to work but use the tag instead.
        if (add_pair(sm,sm2id,fn,first_sm.s) < 0) goto fail;
    }

    ks_free(&first_sm);
    ks_free(&buf);
    ks_free(&id_val);
    ks_free(&sm_val);
    return 0;

 fail:
    ks_free(&first_sm);
    ks_free(&buf);
    ks_free(&id_val);
    ks_free(&sm_val);
    return -1;
}

int bam_smpl_rg2smid(const bam_sample_t *sm, const char *fn, const char *rg, kstring_t *str)
{
    khint_t k;
    khash_t(sm) *rg2smid = (khash_t(sm)*)sm->rg2smid;
    if (rg) {
        str->l = 0;
        kputs(fn, str); kputc('/', str); kputs(rg, str);
        k = kh_get(sm, rg2smid, str->s);
    } else k = kh_get(sm, rg2smid, fn);
    return k == kh_end(rg2smid)? -1 : kh_val(rg2smid, k);
}
