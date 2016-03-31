using System;
using System.Collections.Generic;

namespace SpectraPh.DataStructures
{
    public class LockableDictionary<TK,TV> : Dictionary<TK,TV>
    {
        private object obj = new object();

        public void AddOrUpdate(TK key, Func<TK, TV> addfunc, Func<TK, TV, TV> updatefunc)
        {
            lock (obj)
            {
                if (!ContainsKey(key))
                    Add(key, addfunc(key));
                else
                {
                    this[key] = updatefunc(key, this[key]);
                }
            }
        }
    }
}
