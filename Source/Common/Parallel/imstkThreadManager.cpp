/*
** This file is part of the Interactive Medical Simulation Toolkit (iMSTK)
** iMSTK is distributed under the Apache License, Version 2.0.
** See accompanying NOTICE for details.
*/

#include "imstkThreadManager.h"
#include "imstkLogger.h"

namespace imstk
{
namespace ParallelUtils
{
#ifdef iMSTK_USE_TBB
std::unique_ptr<tbb::global_control> ThreadManager::s_tbbGlobalControl;
#endif

void
ThreadManager::setThreadPoolSize(const size_t nThreads)
{
#ifdef iMSTK_USE_TBB
    LOG_IF(FATAL, (nThreads == 0)) << "Invalid number of threads";
    LOG(INFO) << "Set number of worker threads to " << nThreads;

    if (s_tbbGlobalControl)
    {
        s_tbbGlobalControl.reset();
    }

    s_tbbGlobalControl = std::unique_ptr<tbb::global_control>(
                new tbb::global_control(tbb::global_control::max_allowed_parallelism,
                                        nThreads));
#endif
}

void
ThreadManager::setOptimalParallelism()
{
#ifdef iMSTK_USE_TBB
    setThreadPoolSize(static_cast<size_t>(tbb::info::default_concurrency()));
#endif
}

size_t
ThreadManager::getThreadPoolSize()
{
#ifdef iMSTK_USE_TBB
    return s_tbbGlobalControl->active_value(tbb::global_control::max_allowed_parallelism);
#else
    return 1;
#endif
}
}  // end namespace ParallelUtils
}  // end namespace imstk
