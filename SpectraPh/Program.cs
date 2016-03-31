using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Configuration;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Runtime.Serialization.Formatters.Binary;
using System.Threading;
using System.Threading.Tasks;
using SpectraPh.Actions;
using SpectraPh.DataStructures;

namespace SpectraPh
{
    internal class Program
    {
        public static int CHROMOSOMES;
        public static int PLOIDITY;
        public static int PRINTDETAIL;

        public static SNPMap[] SNPmap;
        public static BlockMatrix Matrix;

        public static long unmappedEnd = 0;
        public static long mappedEnd = 0;

        public static int BINSIZE;

        private static void Main(string[] args)
        {
            var actions = GetAvailableActions();
            if (args.Length == 0)
            {
                PrintUsage(actions);
                Environment.Exit(0);
            }

            if (args.Length >= 1)
            {
                // action was specified
                string actionName = args[0];

                IHaploParserAction action;
                if (!actions.TryGetValue(actionName, out action))
                {
                    Console.WriteLine("Unknown action");
                    PrintUsage(actions);
                    Environment.Exit(1);
                }

                if (args.Length == 1)
                {
                    // only show help on action
                    PrintUsage(action);
                    Environment.Exit(0);
                }

                // if we got here we shoud run the action
                action.Go(args.Skip(1).ToArray());
            }

            Helpers.Write((int)PrintDetail.Important, "Done. Press enter to quit.");
            //Console.ReadLine();
        }

        #region Action Handling
        private static Dictionary<string, IHaploParserAction> GetAvailableActions()
        {
            var results = new Dictionary<string, IHaploParserAction>();
            var type = typeof(IHaploParserAction);
            var actions =
                AppDomain.CurrentDomain.GetAssemblies().SelectMany(s => s.GetTypes()).Where(type.IsAssignableFrom);

            foreach (var action in actions)
            {
                if (!action.IsClass || action.IsAbstract)
                    continue;

                var instance = Activator.CreateInstance(action) as IHaploParserAction;
                results.Add(instance.Command(), instance);
            }

            return results;
        }

        private static void PrintUsage(Dictionary<string, IHaploParserAction> availableActions)
        {
            Console.WriteLine("Usage: ");
            Console.WriteLine("\tHaploSamParser.exe <command> <args>");
            Console.WriteLine("Available commnands:");
            Console.WriteLine("\t" + String.Join("\r\n\t", availableActions.Keys));
        }

        private static void PrintUsage(IHaploParserAction action)
        {
            Console.WriteLine("Wrong usage of the action, please notice it's usage:");
            Console.WriteLine(action.Command() + " " + action.Usage());
        }
        #endregion
    }
}