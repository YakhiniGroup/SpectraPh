using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Reflection;
using System.Text.RegularExpressions;

namespace SpectraPh
{
    public static class Helpers
    {
        readonly static FieldInfo charPosField = typeof(StreamReader).GetField("charPos", BindingFlags.NonPublic | BindingFlags.Instance | BindingFlags.DeclaredOnly);
        readonly static FieldInfo byteLenField = typeof(StreamReader).GetField("byteLen", BindingFlags.NonPublic | BindingFlags.Instance | BindingFlags.DeclaredOnly);
        readonly static FieldInfo charBufferField = typeof(StreamReader).GetField("charBuffer", BindingFlags.NonPublic | BindingFlags.Instance | BindingFlags.DeclaredOnly);
        public readonly static Regex CigarRegex = new Regex(@"([0-9]+[MIDNSHPX=])");

        public static byte TupleId(this Tuple<char, char> tup, char val)
        {
            return (byte) (tup.Item1 == val ? 1 : tup.Item2 == val ? (byte?) 2 : 4);
        }

        public static IEnumerable<string> GetFilesFromPath(string path, string type = "*")
        {
            var attr = File.GetAttributes(path);
            if ((attr & FileAttributes.Directory) == FileAttributes.Directory)
                return Directory.GetFiles(path, type);
            else
                return new List<string>() { path };
        }

        public static void Write(int level, string message)
        {
            if (level >= Program.PRINTDETAIL)
                Console.WriteLine("{0} :: {1}", DateTime.Now, message);
        }

                public static long GetPosition(this StreamReader reader)
        {
            //shift position back from BaseStream.Position by the number of bytes read
            //into internal buffer.
            int byteLen = (int)byteLenField.GetValue(reader);
            var position = reader.BaseStream.Position - byteLen;

            //if we have consumed chars from the buffer we need to calculate how many
            //bytes they represent in the current encoding and add that to the position.
            int charPos = (int)charPosField.GetValue(reader);
            if (charPos > 0)
            {
                var charBuffer = (char[])charBufferField.GetValue(reader);
                var encoding = reader.CurrentEncoding;
                var bytesConsumed = encoding.GetBytes(charBuffer, 0, charPos).Length;
                position += bytesConsumed;
            }

            return position;
        }

        public static void SetPosition(this StreamReader reader, long position)
        {
            reader.DiscardBufferedData();
            reader.BaseStream.Seek(position, SeekOrigin.Begin);
        }

        private static ConcurrentDictionary<string, int> chromNums = new ConcurrentDictionary<string, int>();
        public static int ChromNum(string chromtext)
        {
            int result;
            if (Int32.TryParse(chromtext, out result)) return result;
            if (chromNums.TryGetValue(chromtext, out result))
                return result;
            if (chromtext.ToLowerInvariant().StartsWith("chr"))
                chromtext = chromtext.Substring(3);
            switch (chromtext)
            {
                case "X":
                case "Y":
                    result = chromNums.GetOrAdd(chromtext, t => Program.CHROMOSOMES);
                    break;
                case "M":
                    result = chromNums.GetOrAdd(chromtext, t => -1);
                    break;
                default:
                    result = chromNums.GetOrAdd(chromtext, t => Convert.ToInt32(chromtext));
                    break;
            }
            return result;
        }

        public static IEnumerable<T> GetSubsetById<T>(this IEnumerable<T> lst, int[] idxes)
        {
            int i = 0, j=0;
            foreach (var item in lst)
            {
                if (idxes[j] == i)
                {
                    yield return item;
                    j++;
                    if (j == idxes.Length)
                        break;
                }
                i++;
            }
        }


    }

    public enum PrintDetail { Debug = 0, Detailed = 1, Important = 2 };


}
