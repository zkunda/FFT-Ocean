using UnityEditor;
using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.Windows;
using System.Collections.Generic;
using System.IO;

public class WavesGenerator : MonoBehaviour
{
    public WavesCascade cascade0;
    public WavesCascade cascade1;
    public WavesCascade cascade2;

    // must be a power of 2
    [SerializeField]
    int size = 256;

    [SerializeField]
    WavesSettings wavesSettings;
    [SerializeField]
    bool alwaysRecalculateInitials = false;
    [SerializeField]
    float lengthScale0 = 250;
    [SerializeField]
    float lengthScale1 = 17;
    [SerializeField]
    float lengthScale2 = 5;

    [SerializeField]
    ComputeShader fftShader;
    [SerializeField]
    ComputeShader initialSpectrumShader;
    [SerializeField]
    ComputeShader timeDependentSpectrumShader;
    [SerializeField]
    ComputeShader texturesMergerShader;

    Texture2D gaussianNoise;
    FastFourierTransform fft;
    Texture2D physicsReadback;

    private void Awake()
    {
        Application.targetFrameRate = -1;
        fft = new FastFourierTransform(size, fftShader);
        gaussianNoise = GetNoiseTexture(size);

        cascade0 = new WavesCascade(size, initialSpectrumShader, timeDependentSpectrumShader, texturesMergerShader, fft, gaussianNoise);
        cascade1 = new WavesCascade(size, initialSpectrumShader, timeDependentSpectrumShader, texturesMergerShader, fft, gaussianNoise);
        cascade2 = new WavesCascade(size, initialSpectrumShader, timeDependentSpectrumShader, texturesMergerShader, fft, gaussianNoise);

        InitialiseCascades();

        physicsReadback = new Texture2D(size, size, TextureFormat.RGBAFloat, false);
    }

    int updateCount = 0;
    void InitialiseCascades()
    {
        updateCount = 0;
        float boundary1 = 2 * Mathf.PI / lengthScale1 * 6f;
        float boundary2 = 2 * Mathf.PI / lengthScale2 * 6f;
        cascade0.CalculateInitials(wavesSettings, lengthScale0, 0.0001f, boundary1);
        cascade1.CalculateInitials(wavesSettings, lengthScale1, boundary1, boundary2);
        cascade2.CalculateInitials(wavesSettings, lengthScale2, boundary2, 9999);

        Shader.SetGlobalFloat("LengthScale0", lengthScale0);
        Shader.SetGlobalFloat("LengthScale1", lengthScale1);
        Shader.SetGlobalFloat("LengthScale2", lengthScale2);
    }

    int sliceCount = 64;
    int columnCount = 8;
    int rowCount = 4;
    private void SaveCascadeSequences()
    {
        var current = RenderTexture.active;

        var c0disp = new Texture2D[sliceCount];
        var c1disp = new Texture2D[sliceCount];
        var c2disp = new Texture2D[sliceCount];
        var c0deri = new Texture2D[sliceCount];
        var c1deri = new Texture2D[sliceCount];
        var c2deri = new Texture2D[sliceCount];
        var c0turb = new Texture2D[sliceCount];
        var c1turb = new Texture2D[sliceCount];
        var c2turb = new Texture2D[sliceCount];
        var timeLen = 2.0f;
        float deltaTime = timeLen / sliceCount;
        for (int i = 0; i < sliceCount; ++i)
        {
            var time = i * deltaTime;
            cascade2.CalculateWavesAtTime(time, deltaTime);
            c2deri[i] = CopyRenderTexture(cascade2.Derivatives);
            c2turb[i] = CopyRenderTexture(cascade2.Turbulence);
        }

        timeLen = 8.0f;
        deltaTime = timeLen / sliceCount;
        for (int i = 0; i < sliceCount; ++i)
        {
            var time = i * deltaTime;
            cascade0.CalculateWavesAtTime(time, deltaTime);
            cascade1.CalculateWavesAtTime(time, deltaTime);
            cascade2.CalculateWavesAtTime(time, deltaTime);
            c1deri[i] = CopyRenderTexture(cascade1.Derivatives);
            c1turb[i] = CopyRenderTexture(cascade1.Turbulence);
        }

        timeLen = 16.0f;
        deltaTime = timeLen / sliceCount;
        for (int i = 0; i < sliceCount; ++i)
        {
            var time = i * deltaTime;
            cascade0.CalculateWavesAtTime(time, deltaTime);
            cascade1.CalculateWavesAtTime(time, deltaTime);
            cascade2.CalculateWavesAtTime(time, deltaTime);

            c0disp[i] = CopyRenderTexture(cascade0.Displacement);
            c0deri[i] = CopyRenderTexture(cascade0.Derivatives);
            c0turb[i] = CopyRenderTexture(cascade0.Turbulence);
        }

        object IntersectRayTriangle(Ray ray, Vector3 v0, Vector3 v1, Vector3 v2, bool bidirectional)
        {
            Vector3 ab = v1 - v0;
            Vector3 ac = v2 - v0;

            // Compute triangle normal. Can be precalculated or cached if
            // intersecting multiple segments against the same triangle
            Vector3 n = Vector3.Cross(ab, ac);

            // Compute denominator d. If d <= 0, segment is parallel to or points
            // away from triangle, so exit early
            float d = Vector3.Dot(-ray.direction, n);
            if (d <= 0.0f && (!bidirectional)) return null;

            // Compute intersection t value of pq with plane of triangle. A ray
            // intersects iff 0 <= t. Segment intersects iff 0 <= t <= 1. Delay
            // dividing by d until intersection has been found to pierce triangle
            Vector3 ap = ray.origin - v0;
            float t = Vector3.Dot(ap, n);
            if ((t < 0.0f) && (!bidirectional)) return null;
            //if (t > d) return null; // For segment; exclude this code line for a ray test

            // Compute barycentric coordinate components and test if within bounds
            Vector3 e = Vector3.Cross(-ray.direction, ap);
            float v = Vector3.Dot(ac, e);
            if (v < 0.0f || v > d) return null;

            float w = -Vector3.Dot(ab, e);
            if (w < 0.0f || v + w > d) return null;

            // Segment/ray intersects triangle. Perform delayed division and
            // compute the last barycentric coordinate component
            float ood = 1.0f / d;
            t *= ood;
            v *= ood;
            w *= ood;
            float u = 1.0f - v - w;

            RaycastHit hit = new RaycastHit();

            hit.point = ray.origin + t * ray.direction;
            hit.distance = t;
            hit.barycentricCoordinate = new Vector3(u, v, w);
            hit.normal = Vector3.Normalize(n);

            return hit;
        }

        //Directory.CreateDirectory("Assets/c0displament");
        BlendCycleTextures(c0disp);
        BlendCycleTextures(c0deri);
        BlendCycleTextures(c1deri);
        BlendCycleTextures(c2deri);
        BlendCycleTextures(c0turb);
        BlendCycleTextures(c1turb);
        BlendCycleTextures(c2turb);

        var prefix = "ocean";
        bool saveWaveData = true;
        MergeTurblences(c0turb, c1turb, c2turb);
        SaveDerivatives(prefix + "Deri0", c0deri);
        SaveDerivatives(prefix + "Deri1", c1deri);
        SaveDerivatives(prefix + "Deri2", c2deri);
        var count = sliceCount / 2;

        var wavemapList = new List<byte[]>();

        var tArray = new Texture2D(size * columnCount, size * rowCount, TextureFormat.RGBA32, false);
        var colorData = tArray.GetPixels32(0);
        float tileWidth = lengthScale0;
        var gridSize = tileWidth / size;
        var ray = new Ray();
        ray.direction = Vector3.down;
        var grid = new Vector3[size, size];
        Vector2 offsetX, offsetZ;
        offsetX = offsetZ = Vector2.zero;
        for (int i = 0; i < count; ++i)
        {
            var tex = c0disp[i + count];
            var fData = tex.GetPixels(0);
            for (int j = 0; j < size; ++j)
            {
                for (int k = 0; k < size; ++k)
                {
                    var v = fData[j * size + k];

                    grid[j, k].Set(v.r, v.g, v.b);
                }
            }

            //if (saveWaveData)
            //{
            //    var heightmap = new byte[size * size];
            //    for (int j = 0; j < size; ++j)
            //    {
            //        for (int k = 0; k < size; ++k)
            //        {
            //            ray.origin = new Vector3((k + 0.5f) * gridSize, 100, (j + 0.5f) * gridSize);
            //            object hit = null;
            //            for (int l = j - 5; l <= j + 5; ++l)
            //            {
            //                for (int m = k - 5; m <= k + 5; ++m)
            //                {
            //                    var a = new Vector3(m * gridSize, 0, l * gridSize) + grid[(l + size) % size, (m + size) % size];
            //                    var b = new Vector3((m + 1) * gridSize, 0, l * gridSize) + grid[(l + size) % size, (m + 1 + size) % size];
            //                    var c = new Vector3((m + 1) * gridSize, 0, (l + 1) * gridSize) + grid[(l + 1 + size) % size, (m + 1 + size) % size];
            //                    var d = new Vector3(m * gridSize, 0, (l + 1) * gridSize) + grid[(l + 1 + size) % size, (m + size) % size];
            //                    //if (l < 0)
            //                    //{
            //                    //    a.z -= tileWidth;
            //                    //    b.z -= tileWidth;
            //                    //}
            //                    //else if (l >= size)
            //                    //{
            //                    //    c.z += tileWidth;
            //                    //    d.z += tileWidth;
            //                    //}
            //                    //if (m < 0)
            //                    //{
            //                    //    a.x -= tileWidth;
            //                    //    d.x -= tileWidth;
            //                    //}
            //                    //else if (m >= size)
            //                    //{
            //                    //    b.x += tileWidth;
            //                    //    c.x += tileWidth;
            //                    //}
            //                    float maxX = Mathf.Max(d.x, Mathf.Max(c.x, Mathf.Max(a.x, b.x)));
            //                    float maxZ = Mathf.Max(d.z, Mathf.Max(c.z, Mathf.Max(a.z, b.z)));
            //                    float minX = Mathf.Min(d.x, Mathf.Min(c.x, Mathf.Min(a.x, b.x)));
            //                    float minZ = Mathf.Min(d.z, Mathf.Min(c.z, Mathf.Min(a.z, b.z)));
            //                    if (ray.origin.x > maxX || ray.origin.x < minX || ray.origin.z < minZ || ray.origin.z > maxZ)
            //                        continue;

            //                    hit = IntersectRayTriangle(ray, a, c, b, true);
            //                    if (hit == null)
            //                    {
            //                        hit = IntersectRayTriangle(ray, a, d, c, true);
            //                        if (hit == null)
            //                            continue;
            //                    }

            //                    byte h = (byte)(Mathf.Clamp(((RaycastHit)hit).point.y / 10.0f * 0.5f + 0.5f, 0, 1) * 255.0f);
            //                    heightmap[j * size + k] = h;
            //                    break;
            //                }
            //                if (hit != null)
            //                    break;
            //            }
            //        }
            //    }
            //    waveHeightmapList.Add(heightmap);
            //}
            //SaveHeightmap("h" + i.ToString(), heightmap, size); // debug only

            var column = i % columnCount;
            var row = (rowCount - 1) - i / columnCount;
            var waveData = new byte[size * size]; // heights
            for (int j = 0; j < fData.Length; ++j)
            {
                var v = fData[j];
                offsetX.x = Mathf.Max(v.r, offsetX.x);
                offsetX.y = Mathf.Min(v.r, offsetX.y);
                offsetZ.x = Mathf.Max(v.b, offsetZ.x);
                offsetZ.y = Mathf.Min(v.b, offsetZ.y);
                //var r = Mathf.Clamp(v.r / 10.0f * 0.5f + 0.5f, 0, 1);
                var h = (byte)(255.0f * Mathf.Clamp(v.g / 10.0f * 0.5f + 0.5f, 0, 1));
                //var b = Mathf.Clamp(v.b / 10.0f * 0.5f + 0.5f, 0, 1);
               // var a = Mathf.Clamp(v.a / 10.0f * 0.5f + 0.5f, 0, 1);
                waveData[j] = h;

                var localColumn = j % size;
                var localRow = j / size;
                colorData[(row * size + localRow) * (size * columnCount) + (column * size + localColumn)] = new Color32(h, 0, 0, 1);
            }
            wavemapList.Add(waveData);
        }

        if (saveWaveData)
        {
            using (var writter = new BinaryWriter(System.IO.File.Create("Assets/waveHeightmaps.bytes")))
            {
                int version = 1;
                writter.Write(version);
                writter.Write(count);
                writter.Write(size);
                foreach (var heightmap in wavemapList)
                {
                    writter.Write(heightmap);
                }
            }
        }

        tArray.SetPixels32(colorData);
        tArray.Apply();
        var bytes = tArray.EncodeToPNG();
        System.IO.File.WriteAllBytes("Assets/" + prefix + "Wave.png", bytes);
        RenderTexture.active = current;
        AssetDatabase.SaveAssets();
        AssetDatabase.Refresh();
    }

    void SaveHeightmap(string name, byte[] heightmap, int size)
    {
        var outputTex = new Texture2D(size, size, TextureFormat.ARGB32, false);
        var colorData = new Color32[size * size];
        for (int i = 0; i < size; ++i)
        {
            for (int j = 0; j < size; ++j)
            {
                var v = heightmap[i * size + j];
                colorData[i * size + j] = new Color32(v, v, v, v);
            }
        }
        outputTex.SetPixels32(colorData);
        outputTex.Apply();
        var bytes = outputTex.EncodeToPNG();
        System.IO.File.WriteAllBytes("Assets/" + name + ".png", bytes);
    }

    private void SaveDerivatives(string name, Texture2D texture)
    {
        var outputTex = new Texture2D(texture.width, texture.height, TextureFormat.ARGB32, false);
        var colorData = outputTex.GetPixels(0); 
        var fData = texture.GetPixels(0);
        for (int i = 0; i < texture.width; ++i)
        {
            for (int j = 0; j < texture.height; ++j)
            {
                var ofs = j * texture.width + i;
                var v = fData[ofs];
                var r = (v.r * 0.5f + 0.5f);
                var g = (v.g * 0.5f + 0.5f);
                var b = (v.b * 0.5f + 0.5f);
                var a = (v.a * 0.5f + 0.5f);
                colorData[ofs] = new Color(r, g, b, a);
            }
        }
        outputTex.SetPixels(colorData);
        outputTex.Apply();
        var bytes = outputTex.EncodeToPNG();
        System.IO.File.WriteAllBytes("Assets/" + name + ".png", bytes);
    }

    void MergeTurblences(Texture2D[] turb0, Texture2D[] turb1, Texture2D[] turb2)
    {
        var count = turb0.Length / 2;
        var tArray = new Texture2D(size * columnCount, size * rowCount, TextureFormat.ARGB32, false);
        var colorData = tArray.GetPixels(0);
        for (int i = 0; i < count; ++i)
        {
            var tex0 = turb0[i + count];
            var tex1 = turb1[i + count];
            var tex2 = turb2[i + count];
            var fData0 = tex0.GetPixels(0);
            var fData1 = tex1.GetPixels(0);
            var fData2 = tex2.GetPixels(0);
            var column = i % columnCount;
            var row = (rowCount - 1) - i / columnCount;
            for (int j = 0; j < fData0.Length; ++j)
            {
                var r = fData0[j].r;
                var g = fData1[j].r;
                var b = fData2[j].r;

                var localColumn = j % size;
                var localRow = j / size;
                colorData[(row * size + localRow) * (size * columnCount) + (column * size + localColumn)] = new Color(r, g, b);
            }
        }
        tArray.SetPixels(colorData);
        tArray.Apply();
        var bytes = tArray.EncodeToPNG();
        System.IO.File.WriteAllBytes("Assets/turb.png", bytes);
    }

    private void SaveDerivatives(string name, Texture2D[] textures)
    {
        var count = textures.Length / 2;
        var tArray = new Texture2D(size * columnCount, size * rowCount, TextureFormat.ARGB32, false);
        var colorData = tArray.GetPixels(0);
        for (int i = 0; i < count; ++i)
        {
            var tex = textures[i + count];
            var fData = tex.GetPixels(0);
            var column = i % columnCount;
            var row = (rowCount - 1) - i / columnCount;
            for (int j = 0; j < fData.Length; ++j)
            {
                var v = fData[j];
                var r = (v.r * 0.5f + 0.5f);
                var g = (v.g * 0.5f + 0.5f);
                var b = (v.b * 0.5f + 0.5f);
                var a = (v.a * 0.5f + 0.5f);

                var localColumn = j % size;
                var localRow = j / size;
                colorData[(row * size + localRow) * (size * columnCount) + (column * size + localColumn)] = new Color(r, g, b, a);
            }
        }
        tArray.SetPixels(colorData);
        tArray.Apply();
        var bytes = tArray.EncodeToPNG();
        System.IO.File.WriteAllBytes("Assets/" + name + ".png", bytes);
    }

    private void BlendCycleTextures(Texture2D [] texArray)
    {
        var half = texArray.Length / 2;
        float interval = 1.0f / (1 + half);
        for (int i = 0; i < half; ++i)
        {
            var t0Data = texArray[i].GetPixelData<float>(0);
            var factor = interval * i;
            var tex = texArray[i + half];
            var tData = tex.GetPixelData<float>(0);
            float hMin = float.MaxValue, hMax = float.MinValue;
            for (int j = 0; j < tData.Length; ++j)
            {
                tData[j] = Mathf.Lerp(tData[j], t0Data[j], factor);
                hMin = Mathf.Min(hMin, tData[j]);
                hMax = Mathf.Max(hMax, tData[j]);
            }
            tex.SetPixelData<float>(tData, 0);
            tex.Apply();
        }
    }

    private Texture2D CopyRenderTexture(RenderTexture rt)
    {
        RenderTexture.active = rt;
        Texture2D result = new Texture2D(rt.width, rt.height, TextureFormat.RGBAFloat, false);
        result.ReadPixels(new Rect(0, 0, rt.width, rt.height), 0, 0);
        result.Apply();
        return result;
    }

    private void SaveRenderTexture(string path, RenderTexture rt)
    {
        RenderTexture.active = rt;
        Texture2D result = new Texture2D(rt.width, rt.height, TextureFormat.RGBAFloat, false);
        result.ReadPixels(new Rect(0, 0, rt.width, rt.height), 0, 0);
        result.Apply();
        AssetDatabase.CreateAsset(result, path);
    }

    private void Update()
    {
        if(UnityEngine.Input.GetKey(KeyCode.P))
        {
            SaveCascadeSequences();
        }

        if (alwaysRecalculateInitials)
        {
            InitialiseCascades();
        }

        cascade0.CalculateWavesAtTime(Time.time, Time.deltaTime);
        //if (updateCount == 0)
        {
            cascade1.CalculateWavesAtTime(Time.time, Time.deltaTime);
            cascade2.CalculateWavesAtTime(Time.time, Time.deltaTime);
        }
        ++updateCount;
        RequestReadbacks();
    }
    Texture2D GetNoiseTexture(int size)
    {
        string filename = "GaussianNoiseTexture" + size.ToString() + "x" + size.ToString();
        Texture2D noise = Resources.Load<Texture2D>("GaussianNoiseTextures/" + filename);
        return noise ? noise : GenerateNoiseTexture(size, true);
    }

    Texture2D GenerateNoiseTexture(int size, bool saveIntoAssetFile)
    {
        Texture2D noise = new Texture2D(size, size, TextureFormat.RGFloat, false, true);
        noise.filterMode = FilterMode.Point;
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                noise.SetPixel(i, j, new Vector4(NormalRandom(), NormalRandom()));
            }
        }
        noise.Apply();

#if UNITY_EDITOR
        if (saveIntoAssetFile)
        {
            string filename = "GaussianNoiseTexture" + size.ToString() + "x" + size.ToString();
            string path = "Assets/Resources/GaussianNoiseTextures/";
            AssetDatabase.CreateAsset(noise, path + filename + ".asset");
            Debug.Log("Texture \"" + filename + "\" was created at path \"" + path + "\".");
        }
#endif
        return noise;
    }

    float NormalRandom()
    {
        return Mathf.Cos(2 * Mathf.PI * Random.value) * Mathf.Sqrt(-2 * Mathf.Log(Random.value));
    }

    private void OnDestroy()
    {
        cascade0.Dispose();
        cascade1.Dispose();
        cascade2.Dispose();
    }

    void RequestReadbacks()
    {
        AsyncGPUReadback.Request(cascade0.Displacement, 0, TextureFormat.RGBAFloat, OnCompleteReadback);
    }

    public float GetWaterHeight(Vector3 position)
    {
        Vector3 displacement = GetWaterDisplacement(position);
        displacement = GetWaterDisplacement(position - displacement);
        displacement = GetWaterDisplacement(position - displacement);

        return GetWaterDisplacement(position - displacement).y;
    }

    public Vector3 GetWaterDisplacement(Vector3 position)
    {
        Color c = physicsReadback.GetPixelBilinear(position.x / lengthScale0, position.z / lengthScale0);
        return new Vector3(c.r, c.g, c.b);
    }

    void OnCompleteReadback(AsyncGPUReadbackRequest request) => OnCompleteReadback(request, physicsReadback);

    void OnCompleteReadback(AsyncGPUReadbackRequest request, Texture2D result)
    {
        if (request.hasError)
        {
            Debug.Log("GPU readback error detected.");
            return;
        }
        if (result != null)
        {
            result.LoadRawTextureData(request.GetData<Color>());
            result.Apply();
        }
    }
}
