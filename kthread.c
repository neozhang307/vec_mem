#include <pthread.h>
#include <stdlib.h>
#include <limits.h>

/************
 * kt_for() *
 ************/

struct kt_for_t;

typedef struct {
	struct kt_for_t *t;
	long i;
} ktf_worker_t;

typedef struct kt_for_t {
	int n_threads;
	long n;
	ktf_worker_t *w;
	void (*func)(void*,long,int);
	void *data;
} kt_for_t;

static inline long steal_work(kt_for_t *t)
{
	int i, min_i = -1;
	long k, min = LONG_MAX;
	for (i = 0; i < t->n_threads; ++i)
		if (min > t->w[i].i) min = t->w[i].i, min_i = i;
	k = __sync_fetch_and_add(&t->w[min_i].i, t->n_threads);
	return k >= t->n? -1 : k;
}

static void *ktf_worker(void *data)
{
	ktf_worker_t *w = (ktf_worker_t*)data;
	long i;
	for (;;) {
		i = __sync_fetch_and_add(&w->i, w->t->n_threads);
		if (i >= w->t->n) break;
		w->t->func(w->t->data, i, w - w->t->w);
	}
	while ((i = steal_work(w->t)) >= 0)
		w->t->func(w->t->data, i, w - w->t->w);
	pthread_exit(0);
}

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n)
{
	int i;
	kt_for_t t;
	pthread_t *tid;
	t.func = func, t.data = data, t.n_threads = n_threads, t.n = n;
	t.w = (ktf_worker_t*)alloca(n_threads * sizeof(ktf_worker_t));
	tid = (pthread_t*)alloca(n_threads * sizeof(pthread_t));
	for (i = 0; i < n_threads; ++i)
		t.w[i].t = &t, t.w[i].i = i;
	for (i = 0; i < n_threads; ++i) pthread_create(&tid[i], 0, ktf_worker, &t.w[i]);
	for (i = 0; i < n_threads; ++i) pthread_join(tid[i], 0);
}

/************
 * kt_for_batch() *
 * Lingqi: *
 * every thread would process a batch of data at a time *
 * this design is aim for prosessing PE, as well as other usage *
 * IMPORTANT: *
 * should consider what would happen if data remaining cannot fill up a batch *
 ************/

struct kt_for_t_batch;

typedef struct {
    struct kt_for_t_batch *t;
    long i;
} ktf_worker_t_batch;

typedef struct kt_for_t_batch {
    int n_threads;
    int s_batch;
    long n;
    ktf_worker_t_batch *w;
    void (*func_batch)(void*,long,int,int);
    void *data;
} kt_for_t_batch;

typedef struct batch_pack
{
    void (*func)(void*,long,int);
    void *data;
}batch_pack;

static inline long steal_work_batch(kt_for_t_batch *t, int step)
{
    int i, min_i = -1;
    long k, min = LONG_MAX;
    for (i = 0; i < t->n_threads; ++i)
        if (min > t->w[i].i) min = t->w[i].i, min_i = i;
    k = __sync_fetch_and_add(&t->w[min_i].i, step);
    return k >= t->n? -1 : k;
}

void process_batch(void *data, long start, int batch, int tid)
{
    long i = start;
    batch_pack *pck = (batch_pack*)data;
    for(int j=0; j<batch; j++,i++)
    {
        pck->func(pck->data, i, tid);//compute it in
    }
}

static void *ktf_worker_batch(void *data)
{
    ktf_worker_t_batch *w = (ktf_worker_t_batch*)data;
    long i;
    int batch;
    int step = w->t->n_threads*w->t->s_batch;
    for (;;) {
        i = __sync_fetch_and_add(&w->i, step);
        batch = w->t->n-i < w->t->s_batch? w->t->n-i:w->t->s_batch;
        if(0>=batch) break;
        w->t->func_batch(w->t->data,i, batch,w - w->t->w);
        if(i+step>w->t->n)break;

    }
    while ((i = steal_work_batch(w->t,step)) >= 0)
    {
        batch = w->t->n-i < w->t->s_batch? w->t->n-i:w->t->s_batch;
        batch = 0>batch? 0:batch;
        w->t->func_batch(w->t->data,i,batch,w - w->t->w);
    }
    pthread_exit(0);
}

void kt_for_batch(int n_threads, int batch_size, void (*func)(void*,long,int), void *data, long n)
{
    int i;
    kt_for_t_batch t;
    batch_pack pck;
    pthread_t *tid;
    pck.func = func;
    pck.data = data;
    t.data = &pck, t.n_threads = n_threads, t.n = n;
    
    t.func_batch = process_batch;
    
    t.s_batch = batch_size>1?batch_size:1;
    
    t.w = (ktf_worker_t_batch*)alloca(n_threads * sizeof(ktf_worker_t_batch));
    tid = (pthread_t*)alloca(n_threads * sizeof(pthread_t));
    
    for (i = 0; i < n_threads; ++i)
        t.w[i].t = &t, t.w[i].i = i*batch_size;
    
    for (i = 0; i < n_threads; ++i) pthread_create(&tid[i], 0, ktf_worker_batch, &t.w[i]);
    for (i = 0; i < n_threads; ++i) pthread_join(tid[i], 0);
}

void kt_for_batch2(int n_threads, int batch_size, void (*func)(void*,long,int,int), void *data, long n)
{
    int i;
    kt_for_t_batch t;
    pthread_t *tid;

    t.data = data, t.n_threads = n_threads, t.n = n;
    
    t.func_batch = func;
    
    t.s_batch = batch_size>1?batch_size:1;
    
    t.w = (ktf_worker_t_batch*)alloca(n_threads * sizeof(ktf_worker_t_batch));
    tid = (pthread_t*)alloca(n_threads * sizeof(pthread_t));
    
    for (i = 0; i < n_threads; ++i)
        t.w[i].t = &t, t.w[i].i = i*batch_size;
    
    for (i = 0; i < n_threads; ++i) pthread_create(&tid[i], 0, ktf_worker_batch, &t.w[i]);
    for (i = 0; i < n_threads; ++i) pthread_join(tid[i], 0);
}
/*****************
 * kt_pipeline() *
 *****************/

struct ktp_t;

typedef struct {
	struct ktp_t *pl;
	int64_t index;
	int step;
	void *data;
} ktp_worker_t;

typedef struct ktp_t {
	void *shared;
	void *(*func)(void*, int, void*);
	int64_t index;
	int n_workers, n_steps;
	ktp_worker_t *workers;
	pthread_mutex_t mutex;
	pthread_cond_t cv;
} ktp_t;

static void *ktp_worker(void *data)
{
	ktp_worker_t *w = (ktp_worker_t*)data;
	ktp_t *p = w->pl;
	while (w->step < p->n_steps) {
		// test whether we can kick off the job with this worker
		pthread_mutex_lock(&p->mutex);
		for (;;) {
			int i;
			// test whether another worker is doing the same step
			for (i = 0; i < p->n_workers; ++i) {
				if (w == &p->workers[i]) continue; // ignore itself
				if (p->workers[i].step <= w->step && p->workers[i].index < w->index)
					break;
			}
			if (i == p->n_workers) break; // no workers with smaller indices are doing w->step or the previous steps
			pthread_cond_wait(&p->cv, &p->mutex);
		}
		pthread_mutex_unlock(&p->mutex);

		// working on w->step
		w->data = p->func(p->shared, w->step, w->step? w->data : 0); // for the first step, input is NULL

		// update step and let other workers know
		pthread_mutex_lock(&p->mutex);
		w->step = w->step == p->n_steps - 1 || w->data? (w->step + 1) % p->n_steps : p->n_steps;
		if (w->step == 0) w->index = p->index++;
		pthread_cond_broadcast(&p->cv);
		pthread_mutex_unlock(&p->mutex);
	}
	pthread_exit(0);
}

void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps)
{
	ktp_t aux;
	pthread_t *tid;
	int i;

	if (n_threads < 1) n_threads = 1;
	aux.n_workers = n_threads;
	aux.n_steps = n_steps;
	aux.func = func;
	aux.shared = shared_data;
	aux.index = 0;
	pthread_mutex_init(&aux.mutex, 0);
	pthread_cond_init(&aux.cv, 0);

	aux.workers = (ktp_worker_t*)alloca(n_threads * sizeof(ktp_worker_t));
	for (i = 0; i < n_threads; ++i) {
		ktp_worker_t *w = &aux.workers[i];
		w->step = 0; w->pl = &aux; w->data = 0;
		w->index = aux.index++;
	}

	tid = (pthread_t*)alloca(n_threads * sizeof(pthread_t));
	for (i = 0; i < n_threads; ++i) pthread_create(&tid[i], 0, ktp_worker, &aux.workers[i]);
	for (i = 0; i < n_threads; ++i) pthread_join(tid[i], 0);

	pthread_mutex_destroy(&aux.mutex);
	pthread_cond_destroy(&aux.cv);
}
